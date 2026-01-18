from io import StringIO, TextIOBase
from collections import defaultdict
import re



class SeqEntryReader:
    """
    Simple iterator over SeqEntry file that yields one record at a time.

    
    Usage:
        for se in SeqEntry("seqentryfile.se"):
            ...
    """
    
    def __init__(self, file):
        self.file = file
        self._file = None
        self._should_close = False

    def __iter__(self):
        self._open_file()
        self._activeSeq = None
        self._activeLines = []
        return self

    def __next__(self):
        while True:
            line = self._file.readline()
            if not line:
                # End of file — yield last record if any
                if self._activeSeq:
                    se = SeqEntry.parse(self._activeLines)
                    self._activeSeq = None
                    self._activeLines = []
                    return se
                raise StopIteration

            line = line.rstrip('\n\r')
            line=line.strip()
            seqName=line.split("\t")[0]
            

            # first record; initialize
            if self._activeSeq is None:
                self._activeSeq=seqName
                self._activeLines=[line]
            # new record; safe and start new one
            elif seqName!=self._activeSeq:
                # New record starts
                # Yield previous record
                se = SeqEntry.parse(self._activeLines)
                self._activeSeq=seqName
                self._activeLines = [line]
                return se
            elif line:  # skip empty lines
                self._activeLines.append(line)

    def _open_file(self):
        if self._file is not None:
            return
        if hasattr(self.file, 'readline'):
            # already a file object
            self._file = self.file
            self._should_close = False
        else:
            # assume it's a path
            path = self.file
            if path.endswith(('.gz', '.gzip')):
                import gzip
                self._file = gzip.open(path, 'rt')
            else:
                self._file = open(path, 'r')
            self._should_close = True

    

    def close(self):
        if self._should_close and self._file is not None:
            self._file.close()
            self._file = None

    def __enter__(self):
        self._open_file()
        return self

    def __exit__(self): 
        self.close()


def isssnp(refc,hash,cov,minc,minfreq):
    # refc = A
    # hash =
    if cov==0:
        return False
    if 'A'!= refc:
        ac=hash['A']
        af=float(ac)/float(cov)
        if ac>=minc and af>=minfreq:
            return True
    if 'T' != refc:
        tc=hash['T']
        tf=float(tc)/float(cov)
        if tc>=minc and tf >=minfreq:
            return True
    if 'C'!= refc:
        cc=hash['C']
        cf=float(cc)/float(cov)
        if cc>=minc and cf>=minfreq:
            return True
    if 'G' != refc:
        gc=hash['G']
        gf=float(gc)/float(cov)
        if gc>=minc and gf >=minfreq:
            return True
    return False


class Indel:
    @classmethod
    def parse(cls,e):
        if len(e)!=5:
            raise Exception(f"Cannot parse Indel {e}")
        ref,type,pos,length,count=e[0],e[1],int(e[2]),int(e[3]),float(e[4])
        return Indel(ref,type,pos,length,count)
 
    def __init__(self,ref:str,type:str,pos:int,length:int,count):
        self.ref=ref
        self.type=type # ins or del
        self.pos=pos
        self.length=length
        self.count=count
    
    def __str__(self):
        # ref, type, pos, length, count
        tp=[self.ref, self.type, f"{self.pos}",f"{self.length}",f"{self.count:.2f}"]
        tp="\t".join(tp)
        return tp
    
    def normalize(self, normfactor:float):
        ni=Indel(self.ref,self.type,self.pos,self.length,float(self.count)/normfactor)
        return ni


class SNP:
    @classmethod
    def parse(cls,e):
        if len(e)!=8:
            raise Exception(f"Cannot parse SNP {e}")
        ref,type,pos,refc,ac,tc,cc,gc =e[0],e[1],int(e[2]),e[3],float(e[4]),float(e[5]),float(e[6]),float(e[7])
        return SNP(ref,pos,refc,ac,tc,cc,gc)


    def __init__(self,ref:str,pos:int,refc:str,ac,tc,cc,gc):
        self.ref=ref
        self.pos=pos
        self.refc=refc
        self.ac=ac
        self.tc=tc
        self.gc=gc
        self.cc=cc
    
    def __str__(self):
        # ref, 'snp', pos, refc, ac tc cc gc
        tp=[self.ref, "snp",  f"{self.pos}", self.refc,f"{self.ac:.2f}",f"{self.tc:.2f}",f"{self.cc:.2f}",f"{self.gc:.2f}"] 
        tp="\t".join(tp)
        return tp
    
    def normalize(self,normFactor:float):
        acn=float(self.ac)/normFactor
        tcn=float(self.tc)/normFactor
        ccn=float(self.cc)/normFactor
        gcn=float(self.gc)/normFactor
        ns=SNP(self.ref,self.pos,self.refc,acn,tcn,ccn,gcn)
        return ns
    



class SeqEntry:
    @classmethod

    def getNormalizationFactor(cls, seqEntries: list, minDistance:int):
        totcoverages=[]
        for se in seqEntries:
            # ignore the ends of the entries
            if len(se.cov) <= 2 *minDistance:
                continue
            if minDistance>0:
                tcov=se.cov[minDistance:-minDistance]
                totcoverages.extend(tcov)
            else:
                totcoverages.extend(se.cov)
        if len(totcoverages)==0:
            raise Exception("Unable to normalize; no valid coverage for a single copy gene")
        mean=float(sum(totcoverages))/float(len(totcoverages))
        return mean
        
        

    @classmethod 
    def parse(cls,lines):
        activeName=None
        covar=None
        ambcovar=None
        snplist=[]
        indellist=[]

        for l in lines:
            tmp=l.split("\t")
            sn=tmp[0]
            if activeName is None:
                activeName=sn
            assert sn == activeName
            feature=tmp[1]
            if feature=="ins" or feature == "del":
                indel=Indel.parse(tmp)
                indellist.append(indel)
            elif feature == "snp":
                snp=SNP.parse(tmp)
                snplist.append(snp)
            elif feature =="cov":
                if covar is not None:
                    raise Exception(f"two coverage arrays for sequence {sn}")
                covar = [float(x) for x in tmp[2].split()]
            elif feature =="ambcov":
                if ambcovar is not None:
                    raise Exception(f"two amb coverage arrays for sequence {sn}")
                ambcovar = [float(x) for x in tmp[2].split()]
            else:
                raise Exception(f"Unknown feature {feature}")
        if covar is None:
            raise Exception(f"No coverage for {activeName}")
        if ambcovar is None:
            raise Exception(f"No ambiguous coverage for {activeName}")
        return SeqEntry(activeName,covar,ambcovar,snplist,indellist)
    
    def __init__(self,seqname:str,cov,ambcov,snplist,indellist):
        self.seqname=seqname
        self.cov=cov
        self.ambcov=ambcov
        self.snplist=snplist
        self.indellist=indellist
    
    def __str__(self):
        # cov
        tmp=" ".join([f"{i:.2f}" for i in self.cov])
        tpcov="\t".join([self.seqname,"cov",tmp])
        #ambcov
        tmp=" ".join([f"{i:.2f}"  for i in self.ambcov])
        tpambcov="\t".join([self.seqname,"ambcov",tmp])
        tp=[tpcov,tpambcov]
        for s in self.snplist:
            tp.append(str(s))
        for id in self.indellist:
            tp.append(str(id))
        topr="\n".join(tp)
        return topr
    
    def normalize(self,normfactor:float):
        cov=[float(i)/normfactor for i in self.cov]
        ambcov=[float(i)/normfactor for i in self.ambcov]
        snplist=[]
        for s in self.snplist:
            snplist.append(s.normalize(normfactor))
        indellist=[]
        for i in self.indellist:
            indellist.append(i.normalize(normfactor))
        return SeqEntry(self.seqname,cov,ambcov,snplist,indellist)


class SeqBuilder:
    def __init__(self,seq:str,seqname:str,minmapq:int):
        self.seq=seq
        self.seqlen=len(seq)
        self.seqname=seqname
        self.snpar=[{'A':0,'T':0,'C':0,'G':0} for i in list(seq)]
        self.covar=[0 for i in list(seq)]
        self.ambcovar=[0 for i in list(seq)]
        self.inscol=[]
        self.delcol=[]
        self.minmapq=minmapq

    def __parse_cigar(self,cigar: str):
        """Parse CIGAR string into list of (op, length) tuples."""
        ops = []
        i = 0
        num = ""
        while i < len(cigar):
            if cigar[i].isdigit():
                num += cigar[i]
            else:
                if num:
                    ops.append((cigar[i], int(num)))
                    num = ""
            i += 1
        return ops
    
    def __addcoverage(self,refpos: int,ops,mapq:int):
        rpos=refpos-1
        qpos=0
        ### ref     ATTTAAACCCC---AAAA
        ### que.    ATTT---CCCCTTTAAAA
        ###             D      I
        for op, length in ops:
            if op in ('H', 'S'):  # Hard/soft clip: consumes query only; does not add to coverage
                qpos += length
            elif op == 'I':  # Insertion: consumes query only; does not add to coverage
                qpos += length
            elif op in('D','N'):  # Deletion: consumes reference only; does add to coverage
                for i in range(length):
                    p=rpos+i
                    if p>=self.seqlen:
                        break
                    self.covar[p]+=1
                    if mapq>=self.minmapq:
                        self.ambcovar[p]+=1
                rpos += length
            elif op in ('M', '=', 'X'):  # Match/mismatch: consumes both; adds coverage
                for i in range(length):
                    p=rpos+i
                    if p>=self.seqlen:
                        break
                    self.covar[p]+=1
                    if mapq>=self.minmapq:
                        self.ambcovar[p]+=1
                rpos += length
                qpos += length

    def __addindels(self,refpos:int,ops):
        rpos=refpos-1
        qpos=0
        ### ref     ATTTAAACCCC---AAAA
        ### que.    ATTT---CCCCTTTAAAA
        ###             D      I
        for op, length in ops:
            if op in ('H', 'S'):  # Hard/soft clip: consumes query only; does not add to coverage
                qpos += length
            elif op == 'I':  # Insertion: consumes query only; does not add to coverage
                self.inscol.append((rpos,length))
                qpos += length
            elif op in('D','N'):  # Deletion: consumes reference only; does add to coverage
                self.delcol.append((rpos,length))
                rpos += length
            elif op in ('M', '=', 'X'):  # Match/mismatch: consumes both
                rpos += length
                qpos += length

    def __addSNPs(self,refpos:int,ops,seq:str):
        ### ref     ATTTAAACCCC---AAAA
        ### que.    ATTT---CCCCTTTAAAA
        rpos=refpos-1
        qpos=0
        for op, length in ops:
            if op in ('H', 'S'):  # Hard/soft clip: consumes query only
                qpos += length
            elif op == 'I':  # Insertion: consumes query only
                qpos += length
            elif op in('D','N'):  # Deletion: consumes reference only
                rpos += length
            elif op in ('M', '=', 'X'):  # Match/mismatch: consumes both
                for i in range(length):
                    base = seq[qpos + i]
                    if base in 'ATCG':
                        p=rpos+i
                        if p>=self.seqlen:
                            break
                        self.snpar[p][base]+=1
                rpos += length
                qpos += length
            # Ignore N (skipped reference), P (padding) if present

    
    def addread(self,refpos:int,cigar:str,mapq:int,seq:str):

        ops=self.__parse_cigar(cigar)
        self.__addcoverage(refpos,ops,mapq) # increase coverage; only cigar and mapquality considered
        self.__addindels(refpos,ops)        # add indels; only cigar considered; mapq ignored
        self.__addSNPs(refpos,ops,seq)      # add snps; only cigar considered; mapq ignored
    
    def toSeqEntry(self,mcsnp,mfsnp,mcindel,mfindel):
        snplist=[]
        for i,snp in enumerate(self.snpar):
            refc=self.seq[i]
            cov=self.covar[i]
            if isssnp(refc,snp,cov,mcsnp,mfsnp):
                snpentry=SNP(self.seqname,i+1,refc,snp['A'],snp['T'],snp['C'],snp['G']) # one based snp position
                snplist.append(snpentry)
        
        indellist=[]
        # INSERTIONS
        tmp=defaultdict(int)
        for ins in self.inscol:
            tmp[ins]+=1
        for ins,count in tmp.items():
            pos=ins[0]-1 # position in ins is 1-based but coverage is 0-based
            cov=self.covar[pos]
            insfreq=float(count)/float(cov)
            if count>=mcindel and insfreq>=mfindel:
                id=Indel(self.seqname,"ins",ins[0],ins[1],count)
                indellist.append(id)

        # DELETIONS; kept separate on purpose; in case I want to treat them differentially later
        tmp=defaultdict(int)
        for de in self.delcol:
            tmp[de]+=1
        for de,count in tmp.items():
            pos=de[0]-1 # position in ins is 1-based but coverage is 0-based
            cov=self.covar[pos]
            defreq=float(count)/float(cov)
            if count>=mcindel and defreq>=mfindel:
                id=Indel(self.seqname,"del",de[0],de[1],count)
                indellist.append(id)

        se=SeqEntry(self.seqname,self.covar,self.ambcovar,snplist,indellist)
        return se

        

def load_fasta(fafile):
    entries = {}
    current_header = None
    current_sequence = []
    fh=None
    if  isinstance(fafile, TextIOBase):
        fh=fafile
    else:
        fh=open(fafile,'r')
    
    for line in fh:
        line = line.rstrip()  # remove trailing \n
        
        if line.startswith('>'):
            # Save previous entry if exists
            if current_header is not None:
                seq = ''.join(current_sequence)
                entries[current_header]=seq
            
            # Start new entry; get rid of >
            current_header = line[1:]
            # split and get rid of anything after whitespace
            if re.search(r'\s', current_header):
                current_header=re.split(r'\s+', current_header)[0]
            current_sequence = []
        elif line and current_header is not None:
            # Add sequence line (skip empty lines)
            current_sequence.append(line)
    
    # Don't forget the last entry!
    if current_header is not None:
        seq = ''.join(current_sequence)
        entries[current_header]=seq
    fh.close()
    
    return entries

def test_computeNormalization():
    ses=[SeqEntry("t",[10,]*10,[],[],[]),SeqEntry("t",[2,]*10,[],[],[])]
    nf=SeqEntry.getNormalizationFactor(ses,0)
    assert nf==6, "test1"

    ses=[SeqEntry("t",[1,10,10,10,10,10,1],[],[],[]),SeqEntry("t",[1,2,2,2,2,2,1],[],[],[])]
    nf=SeqEntry.getNormalizationFactor(ses,0)
    assert nf<5, "test2"
    nf=SeqEntry.getNormalizationFactor(ses,1)
    assert nf==6, "test3"

    print("Quick test computation of normalization factor passed ✓")


def test_normalize():
    s=SNP("chr1",1,"A",5,6,7,1)
    sn=s.normalize(2.0)
    assert sn.ref=="chr1"
    assert sn.pos== 1
    assert sn.refc=="A"
    assert sn.ac==2.5
    assert sn.tc==3
    assert sn.cc==3.5
    assert sn.gc==0.5

    print("Quick test of SNP normalization PASSED ✓")
    id=Indel("chr2", "ins",5,2,11)
    idn=id.normalize(2.0)
    assert idn.ref == "chr2"
    assert idn.type=="ins"
    assert idn.pos==5
    assert idn.count==5.5
    assert idn.length==2
    print("Quick test of Insertion normalization PASSED ✓")
    deli=Indel("chr3", "del",5,2,20)
    de=deli.normalize(5.0)
    assert de.ref == "chr3"
    assert de.type=="del"
    assert de.pos==5
    assert de.count==4
    assert de.length==2
    print("Quick test of Deletion normalization PASSED ✓")

    id=Indel("chr2", "ins",5,2,11)
    deli=Indel("chr3", "del",5,2,20)
    s=SNP("chr1",1,"A",5,6,7,1)
    se=SeqEntry("te1",[5,6,6,4,2],[2,3,4,6,1],[s],[id,deli])
    sn=se.normalize(2)
    assert sn.cov[0]==2.5
    assert sn.cov[1]==3
    assert sn.cov[4]==1
    assert sn.ambcov[0]==1
    assert sn.ambcov[1]==1.5
    assert sn.ambcov[4]==0.5
    assert sn.ambcov[3]==3
    assert sn.snplist[0].ac==2.5
    assert sn.indellist[0].count==5.5
    assert sn.indellist[1].count==10
    print("Quick test of SeqEntry normalization PASSED ✓")




def test_Seq_Builder_add():
    sb=SeqBuilder("AAATTTCCCGGG","hans",5)
    sb.addread(1,"3M",4,"AAA")
    sb.addread(1,"3M",5,"TTT")

    assert sb.covar[0]==2
    assert sb.ambcovar[0]==1
    assert sb.covar[1]==2
    assert sb.ambcovar[1]==1
    assert sb.covar[2]==2
    assert sb.ambcovar[2]==1
    assert sb.covar[3]==0
    assert sb.ambcovar[3]==0
    assert sb.snpar[0]['A']==1
    assert sb.snpar[0]['T']==1
    

    # AAATTT---CCCGGG
    # 123456---789012
    #    TTTAAACCC
    sb.addread(4,"3=3I3X",5,"TTTAAACCC")
    assert sb.covar[3]==1
    assert sb.covar[4]==1
    assert sb.covar[5]==1
    assert sb.covar[6]==1
    assert sb.covar[7]==1
    assert sb.covar[8]==1
    assert sb.covar[9]==0
    assert sb.snpar[3]['T']==1
    assert sb.snpar[6]['A']==0
    assert sb.snpar[6]['C']==1
    assert sb.inscol[0]==(6,3),  f"got {sb.inscol[0]}"


    # AAATTTCCCGGG
    # 123456789012
    #    TTT---AAA
    sb.addread(4,"3M3D3M",5,"TTTAAA")
    assert sb.covar[3]==2
    assert sb.covar[4]==2
    assert sb.covar[5]==2
    assert sb.covar[6]==2
    assert sb.covar[7]==2
    assert sb.covar[8]==2
    assert sb.covar[9]==1
    assert sb.covar[10]==1
    assert sb.covar[11]==1
    assert sb.delcol[0]==(6,3), f"got {sb.delcol[0]}"

    sb.addread(12,"3M",5,"TTT")


    print("Quick test of SeqBuilder add PASSED ✓")


def test_Seq_Builder_init():
    sb=SeqBuilder("AAATTTCCCGGG","hans",5)
    assert sb.seq == "AAATTTCCCGGG",           f"sequence"
    assert sb.seqname == "hans",            f"seqname"
    assert sb.minmapq ==5,                  f"minmapq"
    assert len(sb.covar) == 12,             f"length of covar"
    assert len(sb.ambcovar) == 12,           f"length of ambcovar"
    assert len(sb.snpar) == 12,           f"length of ambcovar"
    assert len(sb.inscol) == 0,           f"length of ambcovar"
    assert len(sb.delcol) == 0,           f"length of ambcovar"
    print("Quick test of SeqBuilder __init__ PASSED ✓")

def test_fasta_loader():
    from io import StringIO
    
    test_content = """>seq1 some description
ACGTACGT
GCTA
>seq2
NNNNNNNNNN
>seq3 empty sequence

>seq4
ATGCATGCATGC
"""

    result = load_fasta(StringIO(test_content))
    
    expected = {
        "seq1": "ACGTACGTGCTA",
        "seq2": "NNNNNNNNNN",
        "seq3": "",
        "seq4": "ATGCATGCATGC"
    }
    
    assert len(result) == 4,           f"Expected 4 sequences, got {len(result)}"
    assert "seq3" in result,           "Missing empty sequence entry"
    assert result["seq3"] == "",       "Empty sequence should be empty string"
    assert result == expected,         "Dictionary content doesn't match expected"
    
    print("Quick test of fasta_loader PASSED ✓")

if __name__ == "__main__":
    test_fasta_loader()
    test_Seq_Builder_init()
    test_Seq_Builder_add()
    test_normalize()
    test_computeNormalization()