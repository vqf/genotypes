# Genotypes
loc_genotype and get_genotype_bam extract and count reads spanning a given genomic region. The main aim of the scripts is to study subclonalities from high-depth sequencing experiments. For instance, this may serve to characterize *in vivo* CRISPR experiments. Small insertions/deletions are highlighted (deletions as -, insertions in parentheses). Reads with hard clipping, reads with soft-clipping spanning the target region and reads not spanning the whole target region are ignored. Both scripts use Samtools from the command line. If the Samtools executable cannot be accessed as `samtools`, the path to the executable can be provided as an argument.

# loc_genotype
This script extracts and classifies the parts of reads that span the target region.

## Input
>perl loc_genotype.pl bam_file chr from to [max_number_of_reads=0] [samtools_executable=samtools]

The input bam_file must be indexed.

## Output

By default, the ouput is written to stdout. Example:

```
INDEL_DEL_GCAC-------------------------------ACCCCCGAGCCGCTGCAGTGGGAAC  1       -31     1
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGT-CCCCGAGCCGCTGCAGTGGGAAC  6       -1      1
GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCAAGCCGCTGCAGTGGGAAC    4       0       0
GCATGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    5       0       0
GCACGGTGCGCGAGCGCAGGTTGTACTCGGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
TCACGGTTCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTAGGAAC    3       0       0
INDEL_DEL_GCACGGTGCGT-----------------------------------GCTGCAGTGGGAAC  4       -35     1
GCACGGTACGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    4       0       0
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTACCTAGCGGGGT-CCCCGAGCCGCTGCAGGGGGAAC  1       -1      1
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGTGT--CCCGAGCCGCTGCAGTGGGAAC  1       -2      1
GCACGGTGCGTGAGCGCAGGGTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCTCCCGAGCCGCTGCAGTGGGAAC    1       0       0
GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCATTGGGAAC    2       0       0
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTACTCAGC-GGGTCCCCCGAGCCGCTGCAGTGGGAAC  37      -1      1
GTACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGCGCCGCTGCAGTGGGAAC    1       0       0
GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGGAC    7       0       0
GCACGGTGCGTGTGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAT    1       0       0
GCACGGTTCGTGAGCGCAGGTTTTACCTAGCCGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTA-------------CCCGAGCCGCTGCAGGGGGAAC  2       -13     1
INDEL_DEL_GCACGATGCGTGAGCGCAGGTTGTACTCAG----------CGAGCCGCTGCAGTGGGAAC  1       -10     1
GCACGGTGCGTGCGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTACTCA-----------------GCTGCAGTGGGAAC  3       -17     1
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTACTCAGC--GGTCCCCGGAGCCGCTGCAGTGGGAAC  1       -2      1
GCACGGTGCGTGAGCGCAGGTCGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    8       0       0
GCACGGTGCGTGGGCGCAGGTTTTACTCGGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGG------------CGCTGCAGTGGGAAC  3       -12     0
GCACGGCGCGTGAGCGCAGGTTGTACTCGGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
GCACGGTGCGTGAGCGCAGGTTGTACGCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
GCACGGTGTGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    3       0       0
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTA--------------CCGAGCCGCTGCAGTGGGAAC  1       -14     1
GCACGGTGCGTTAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    1       0       0
GCACGGTGCGTGAGCGCAGGTTGTACTCGGCGGGGTCCCCCGAGCCGCTGCAGTGGGAAC    9       0       0
INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTACTCAG-------CCCCGAGCCGCTGCAGTGGGAAC  1       -7      1
INDEL_DEL_--------------------------------------------------------GAAC  3       -56     1
INDEL_DEL_GCACGGTGCGTGAGCGCA---------------GGTCCCCCGAGTCGCTGCAGTGGGAAC  1       -15     0
GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCACTGCAGTGGGAAC    2       0       0
INDEL_INS_GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGG(TGCGGGGCAGCATCAACAGATTCAAGACCAGCGACTACGTGAAAGAAGCCAAACAGCTGCTGAA)GGTCCCCCGAGCCGCTGCAGTGGGAAC        1       64      1
```
Reads featuring small insertions/deletions and soft-clipping have possibly combined tags at the beginning (`INDEL_INS_`, `INDEL_DEL_` and `S_`). The columns are `sequence`, `number_of_reads`, `delta_frame` and `frameshift`. `delta_frame` keeps track of insertions and deletions in the region, and `frameshift` takes a 0 if those changes do not affect the reading frame and 1 otherwise.

# get_genotype_bam

This script retrieves all the reads with a given sequence at the target region (as provided by loc_genotype) and writes them as a bam file.
## Input
>perl get_genotype_bam.pl bam_file chr from to type_file [samtools_executable=samtools]

- `bam_file` must be indexed. 
- type_file is the name of a file containing the expected sequences in the target region (`chr`:`from`-`to`) as provided by `loc_genotype`. Example,
```
  INDEL_DEL_GCACGGTGCGTGAGCGCAGGTTGTACTCAGC-GGGTCCCCCGAGCCGCTGCAGTGGGAAC
  GCACGGTGCGTGAGCGCAGGTTGTACTCAGCGGGGTCCCCCGAGCCGCTGCAGTGGGGAC
```
## Output
This `type_file` would create two bam files: bam_0.bam with 37 reads (containing the first target sequence) and bam_1.bam with 2 reads.