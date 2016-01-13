python bwa_parallel.py -P tmp -I . -R /prj/gf-rhizo/data/reference/Bm_WB42_1.2/concat/Bm_WB42_1.2.a_concat.fasta -t 4 -e error -o out -n 50
If the .fastq files are apropriately named, it is sufficient to just call bwa_parallel.py, and tell it the location of the fatsq files -I, the location of the reference to map to -R, and the mapping parameters -t,-n.
