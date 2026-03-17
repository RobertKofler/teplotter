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
        tmp=PlotableFormater.prepareForPrint(se,args.sampleid,tomask,args.ymax)
        tr=[format_col(i) for i in tmp] # final formatting step; padding
        tp="\n".join(tr)

        if args.outputdir is not None:
            filename=se.seqname
            filename=filename.replace("/","_")
            filename=filename.replace(" ","_")
            full_path = os.path.join(args.outputdir, filename)
            with open(full_path, "w") as f:
                f.write(tp)
        else:
            writer.write(tp)







