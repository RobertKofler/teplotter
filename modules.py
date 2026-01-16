from io import StringIO, TextIOBase
import re

class SeqBuilder:
    def __init__(self,seq:str,seqname:str,minmapq:int):
        self.seq=seq
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
                    self.covar[rpos+i]+=1
                if mapq>=self.minmapq:
                    for i in range(length):
                        self.ambcovar[rpos+i]+=1
                rpos += length
            elif op in ('M', '=', 'X'):  # Match/mismatch: consumes both; adds coverage
                for i in range(length):
                    self.covar[rpos+i]+=1
                if mapq>=self.minmapq:
                    for i in range(length):
                        self.ambcovar[rpos+i]+=1
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
                self.inscol.append([rpos+1,length])
                qpos += length
            elif op in('D','N'):  # Deletion: consumes reference only; does add to coverage
                self.delcol.append([rpos+1,length])
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
            elif op == 'D':  # Deletion: consumes reference only
                rpos += length
            elif op in ('M', '=', 'X'):  # Match/mismatch: consumes both
                for i in range(length):
                    base = seq[qpos + i]
                    if base in 'ATCG':
                        self.snpar[rpos+i][base]+=1
                rpos += length
                qpos += length
            # Ignore N (skipped reference), P (padding) if present

    
    def addread(self,refpos:int,cigar:str,mapq:int,seq:str):

        ops=self.__parse_cigar(cigar)
        self.__addcoverage(refpos,ops,mapq) # increase coverage; only cigar and mapquality considered
        self.__addindels(refpos,ops)        # add indels; only cigar considered; mapq ignored
        self.__addSNPs(refpos,ops,seq)      # add snps; only cigar considered; mapq ignored

        

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