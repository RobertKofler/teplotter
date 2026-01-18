#!/usr/bin/env python
import argparse
import modules
import os




parser = argparse.ArgumentParser(description="""           
normalizes the coverage for seqentries
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
parser.add_argument('--seq-entries', type=str, default=None,dest="seqentry", required=True, help="A seq entry file")
parser.add_argument("--scg-ids", type=str, required=True, dest="seqids", default=None, help="IDs of the entries that should be plotted; separated by comma; can also be 'ALL'")
parser.add_argument("--sample-id", type=str, required=False, dest="sampleid", default="x", help="the ID of current sample")
parser.add_argument("--output-dir", type=str, required=False, dest="outputdir", default=None, help="the output directory")

args = parser.parse_args()
seqset=None
if "," in args.seqids:
    seqset=set(args.seqids.split(","))
else:
    seqset=set([args.seqids])
printall=False
if args.seqids.lower() == "all":
    True



args = parser.parse_args()
for se in modules.SeqEntryReader(args.seqentry):
    if printall or se.seqname in seqset:
        tp=prepareForPrint(se,args.sampleid)
        if args.outputdir is not None:
            full_path = os.path.join(args.outputdir, se.seqname)
            with open(full_path, "w") as f:
                f.write(tp)
        else:
            print(tp)



def prepareForPrint(se:modules.SeqEntry, sampleid:str):
    lines=[]
    for i,c in enumerate(se.cov):
        # seqname, sampleid, cov, pos, count
        tmp=[se.seqname,sampleid,"cov",str(i+1),str(c)]
        lines.append("\t".join(tmp))

    for i,ac in enumerate(se.ambcov):
        # seqname, sampleid, ambcov, pos, count
        tmp=[se.seqname,sampleid,"ambcov",str(i+1),str(ac)]
        lines.append("\t".join(tmp))

    for s in se.snplist:
        # seqname, sampleid, snp, pos, refc, ac, tc, cc, gc
        # SNP(ref,pos,refc,ac,tc,cc,gc)
        a={"A":s.ac,"T":s.tc,"C":s.cc,"G":s.gc}
        a[s.refc]=0 # do not visualize the reference allele for a SNP
        tmp=[se.seqname,sampleid,"snp",str(s.pos), s.refc,str(a.ac) , str(a.tc) , str(a.cc) , str(a.gc)]
        lines.append("\t".join(tmp))
    
    for i in se.indellist:
        if i.type=="del":
            # seqname, sampleid, del, pos, length, count
            tmp=[se.seqname,sampleid,"del",str(i.pos),str(i.length),str(i.count)]
            lines.append("\t".join(tmp))
            # ref:str,type:str,pos:int,length:int,count

        elif i.type=="ins":
            # seqname, sampleid, ins, startpos, endpos, startcov,endcov, count
                # AAATTTCCCGGG
                # 123456789012
                #    TTT---AAA
                # pos = 6 and len = 3
                # bow from 6 to 10 (actual 0-based coverages are 5 and 9)
            startpos=i.pos
            endpos=startpos+i.length+1
            startcov=se.cov[startpos-1]
            endcov=se.cov[endpos-1]
            tmp=[se.seqname,sampleid,"ins",str(i.startpos),str(i.endpos),str(startcov),str(endcov),str(i.count)]
            lines.append("\t".join(tmp))

        else:
            raise Exception(f"invalid type{i.type}")
        
    tr="\n".join(lines)
    return tr



