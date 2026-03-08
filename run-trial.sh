samtools view ../../dev/data/Dmel01_scg_te_library.sorted.bam|python bam2so.py --sam - --fasta ../../dev/data/Dmel_te_scg.fasta > ../../dev/data/Dmel01.summary
python normalize-so.py --so ../../dev/data/Dmel01.summary --end-distance 100 --exclude-quantile 25 > ../../dev/data/Dmel01.normalized
