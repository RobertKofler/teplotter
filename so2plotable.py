#!/usr/bin/env python
import argparse
import os
from collections import defaultdict
from modules import PlotableFormater,Writer,SeqEntryReader, load_bed






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





def prepareForPrint(se:SeqEntry, sampleid:str,tomask,ymax):
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
    covt=prepareCoveragForPrint(s ,cov,sampleid,"cov")
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
        

    return lines


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
writer=Writer(args.outfile)

tomask = load_bed(args.maskbed) # 0-based; bed is 0-based


seqset=None
if "," in args.seqids:
    seqset=set(args.seqids.split(","))
else:
    seqset=set([args.seqids])
printall=False
if args.seqids.lower() == "all":
    printall=True


for se in SeqEntryReader(args.so):
    if printall or se.seqname in seqset:
        tmp=prepareForPrint(se,args.sampleid,tomask,args.ymax)
        tr=[format_col(i) for i in tmp] # final formatting step; padding
        tp="\n".join(tp)

        if args.outputdir is not None:
            filename=se.seqname
            filename=filename.replace("/","_")
            filename=filename.replace(" ","_")
            full_path = os.path.join(args.outputdir, filename)
            with open(full_path, "w") as f:
                f.write(tp)
        else:
            writer.write(tp)







