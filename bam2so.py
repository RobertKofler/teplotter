#!/usr/bin/env python
import argparse
import modules




parser = argparse.ArgumentParser(description="""           
summarize coverage for diverse features
""",formatter_class=argparse.RawDescriptionHelpFormatter,
epilog="""
Authors
-------
    Robert Kofler
""")
parser.add_argument('--sam', type=argparse.FileType('r'), default=None,dest="sam", required=True, help="A sam file")
parser.add_argument("--fasta", type=str, required=True, dest="fasta", default=None, help="the fasta file to which reads were mapped")
parser.add_argument("--mapqth", type=int, required=False, dest="mapqth", default=5, help="mapping quality threshold; below ambiguous")
parser.add_argument("--mc-snp", type=int, required=False, dest="mcsnp", default=5, help="minimum count of SNPs")
parser.add_argument("--mf-snp", type=float, required=False, dest="mfsnp", default=0.1, help="minimum frequency of SNPs")
parser.add_argument("--mc-indel", type=int, required=False, dest="mcindel", default=3, help="minimum count of indels")
parser.add_argument("--mf-indel", type=float, required=False, dest="mfindel", default=0.01, help="minimum frequency of indels")
parser.add_argument("--output-file", type=str, required=False, dest="outfile", default=None, help="output file in so format; if none is provided output will be screen")

args = parser.parse_args()
writer=modules.Writer(args.outfile)

# load fasta from file into dict
fastalib=modules.load_fasta(args.fasta)

#
activeBuilder=None

for line in args.sam:
    line = line.strip()
    if not line or line.startswith('@'):
        continue  # Skip header
    
    fields = line.split('\t')
    if len(fields) < 11:
        continue  # Malformed
    flag = int(fields[1])
    if flag & 0x4:  # Unmapped
        continue
    if flag & 0x100:  # Secondary alignment
        continue
    if flag & 0x800:  # Supplementary
        continue

    # initialize fields
    ref = fields[2]    
    pos = int(fields[3])  # 1-based start position
    mapq = int(fields[4])
    cigar = fields[5]
    seq = fields[9].upper()  # Query sequence
    if cigar == '*':
            continue
    
    # start of script
    if activeBuilder is None:
         refseq=fastalib[ref]
         activeBuilder=modules.SeqBuilder(refseq,ref,args.mapqth)
    
    # new refseq
    if ref != activeBuilder.seqname:
         sbe=activeBuilder.toSeqEntry(args.mcsnp,args.mfsnp,args.mcindel,args.mfindel)
         writer.write(str(sbe))
         # print entry
         refseq=fastalib[ref]
         activeBuilder=modules.SeqBuilder(refseq,ref,args.mapqth)
    
    # in any case add the read
    activeBuilder.addread(pos,cigar,mapq,seq)

# process the last one as well
sbe=activeBuilder.toSeqEntry(args.mcsnp,args.mfsnp,args.mcindel,args.mfindel)
writer.write(str(sbe))



 