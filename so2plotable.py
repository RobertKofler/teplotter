#!/usr/bin/env python
import argparse
import modules
import os

padto=9
def format_col(topr:list):
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


def prepareForPrint(se:modules.SeqEntry, sampleid:str):
    lines=[]
    covt=prepareCoveragForPrint(se.cov,sampleid,"cov")
    ambcovt=prepareCoveragForPrint(se.ambcov,sampleid,"ambcov")
    lines.extend(covt)
    lines.extend(ambcovt)



    for s in se.snplist:
        # seqname, sampleid, snp, pos, refc, ac, tc, cc, gc
        # SNP(ref,pos,refc,ac,tc,cc,gc)
        a={"A":s.ac,"T":s.tc,"C":s.cc,"G":s.gc}
        a[s.refc]=0 # do not visualize the reference allele for a SNP
        tmp=[se.seqname,sampleid,"snp",str(s.pos), s.refc,str(a["A"]) , str(a["T"]) , str(a["C"]) , str(a["G"])]
        lines.append(format_col(tmp))
    
    for i in se.indellist:
        if i.type=="del":
            # seqname, sampleid, del, pos, length, count
            tmp=[se.seqname,sampleid,"del",str(i.pos),str(i.length),str(i.count)]
            lines.append(format_col(tmp))
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
            tmp=[se.seqname,sampleid,"ins",str(startpos),str(endpos),str(startcov),str(endcov),str(i.count)]
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



args = parser.parse_args()

if args.outfile is not None and args.outputdir is not None:
    raise Exception("invalid parameters; either provide output-dir or output-file; not both")
# initialize writer
writer=modules.Writer(args.outfile)

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
        tp=prepareForPrint(se,args.sampleid)
        if args.outputdir is not None:
            filename=se.seqname
            filename=filename.replace("/","_")
            filename=filename.replace(" ","_")
            full_path = os.path.join(args.outputdir, filename)
            with open(full_path, "w") as f:
                f.write(tp)
        else:
            writer.write(tp)







