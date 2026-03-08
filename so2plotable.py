#!/usr/bin/env python
import argparse
import modules
import os
from collections import defaultdict

def readbed(bed_path: str):
    """
    Reads a BED file and returns a dictionary of positions that are present in the file.
    bed is 0-based coordinates
    """

    # Using defaultdict(dict) for clean nested structure
    result = defaultdict(lambda: defaultdict(bool))
    if bed_path is None:
        return result
    
    with open(bed_path, 'rt') as f:   
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#') or line.startswith('track ') or line.startswith('browser '):
                continue
                
            fields = line.split('\t')
            assert len(fields)>=3
            chrom = fields[0]
            start = int(fields[1])
            end   = int(fields[2])

            # BED is [start, end) → we include all positions from start inclusive to end inclusive
            for pos in range(start, end+1):
                result[chrom][pos] = True
    return result




padto=9
def format_col(topr:list):
    """
    format a single line and pad empty entries; returns string
    
    :param topr: list of entries 
    :type topr: list
    """
    if padto==0:
        return "\t".join(topr)
    tf=[]
    for i in range(padto):
        if(i< len(topr)):
            tf.append(str(topr[i]))
        else:
            tf.append("")
    return "\t".join(tf)

def prepareCoveragForPrint(set:list, sampleid:str,covtype:str):
    """
    print the coverage a list of coverages
    
    :param set: the coverages to be printed
    :type set: list
    :param sampleid: the sample id
    :type sampleid: str
    :param covtype: the type of coverage; eg cov, ambcov, mcov 
    :type covtype: str
    """
    tmp=[]
    for i,c in enumerate(set):
        # seqname, sampleid, cov, pos, count
        t=[se.seqname,sampleid,covtype,str(i+1),str(c)]
        tmp.append(t)
       
    first,last=tmp[0],tmp[-1]
    newfirst=[first[0],first[1],first[2],first[3],"0.0"]
    newlast=[last[0],last[1],last[2],last[3],"0.0"]
    tmp.insert(0,newfirst)
    tmp.append(newlast)

    topr=[]
    for i in tmp:
        topr.append(format_col(i))
    return topr


def prepareForPrint(se:modules.SeqEntry, sampleid:str,tomask,ymax):
    # get local masking
    localmask=tomask[se.seqname] # bed is 0-based
    # coverages and mask according to user specifications
    cov=se.cov
    ambcov=se.ambcov
    mcov=[0]*len(cov)
 
    for i in range(0,len(cov)):
        c=cov[i]
         # mask coverage if either in localmaks or coverage exceeds ymax
        if (i in localmask):
            ambcov[i]=0
            mcov[i]=cov[i]
            cov[i]=0
        elif ymax is not None and c>ymax:
            ambcov[i]=0
            mcov[i]=ymax
            cov[i]=0
            localmask[i]=True

    lines=[]
    covt=prepareCoveragForPrint(cov,sampleid,"cov")
    ambcovt=prepareCoveragForPrint(ambcov,sampleid,"ambcov")
    mcovt=prepareCoveragForPrint(mcov,sampleid,"mcov")
    lines.extend(covt)
    lines.extend(ambcovt)
    lines.extend(mcovt)



    for s in se.snplist:
        if s.pos in localmask:
            continue
        # seqname, sampleid, snp, pos, refc, ac, tc, cc, gc
        # SNP(ref,pos,refc,ac,tc,cc,gc)
        a={"A":s.ac,"T":s.tc,"C":s.cc,"G":s.gc}
        for base,count in a.items():
            if count ==0 or base==s.refc:
                continue
            tmp=[se.seqname,sampleid,"snp",str(s.pos), s.refc,base,str(count)]
            lines.append(format_col(tmp))
    
    for i in se.indellist:


        if i.type=="ins":
            if i.pos in localmask:
                continue
            # seqname, sampleid, del, pos, length, count
            tmp=[se.seqname,sampleid,"ins",str(i.pos),str(i.length),str(i.count)]
            lines.append(format_col(tmp))
            # ref:str,type:str,pos:int,length:int,count

        elif i.type=="del":
            # seqname, sampleid, ins, startpos, endpos, startcov,endcov, count
                # AAATTTCCCGGG
                # 123456789012
                #    TTT---AAA
                # pos = 6 and len = 3
                # bow from 6 to 10 (actual 0-based coverages are 5 and 9)
            startpos=i.pos
            endpos=startpos+i.length+1
            if startpos in localmask or endpos in localmask:
                continue
            startcov=se.cov[startpos-1]
            endcov=se.cov[endpos-1]
            tmp=[se.seqname,sampleid,"del",str(startpos),str(endpos),str(startcov),str(endcov),str(i.count)]
            lines.append(format_col(tmp))

        else:
            raise Exception(f"invalid type{i.type}")
        
    tr="\n".join(lines)
    return tr


parser = argparse.ArgumentParser(description="""           
normalizes the coverage for seqentries
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
parser.add_argument('--so', type=str, default=None,dest="so", required=True, help="A sequence overview (so) file")
parser.add_argument("--seq-ids", type=str, required=True, dest="seqids", default=None, help="IDs of the entries that should be plotted; separated by comma; can also be 'ALL'")
parser.add_argument("--sample-id", type=str, required=False, dest="sampleid", default="x", help="the ID of current sample")
parser.add_argument("--output-dir", type=str, required=False, dest="outputdir", default=None, help="the output directory; a plotable will be written for each fasta entry")
parser.add_argument("--output-file", type=str, required=False, dest="outfile", default=None, help="output file in plotable format;")
parser.add_argument("--mask-bed", type=str, required=False, dest="maskbed", default=None, help="a bed file for masking; regions in the file will be masked")
parser.add_argument("--mask-ymax", type=int, required=False, dest="ymax", default=None, help="mask regions with coverages exceeding the given value")

args = parser.parse_args()

if args.outfile is not None and args.outputdir is not None:
    raise Exception("invalid parameters; either provide output-dir or output-file; not both")
# initialize writer
writer=modules.Writer(args.outfile)

tomask = readbed(args.maskbed) # 0-based; bed is 0-based


seqset=None
if "," in args.seqids:
    seqset=set(args.seqids.split(","))
else:
    seqset=set([args.seqids])
printall=False
if args.seqids.lower() == "all":
    printall=True


for se in modules.SeqEntryReader(args.so):
    if printall or se.seqname in seqset:
        tp=prepareForPrint(se,args.sampleid,tomask,args.ymax)
        if args.outputdir is not None:
            filename=se.seqname
            filename=filename.replace("/","_")
            filename=filename.replace(" ","_")
            full_path = os.path.join(args.outputdir, filename)
            with open(full_path, "w") as f:
                f.write(tp)
        else:
            writer.write(tp)







