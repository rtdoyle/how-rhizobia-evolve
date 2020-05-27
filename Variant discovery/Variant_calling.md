Whole genome sequencing, align reads, call variants
================
Rebecca Batstone
2020-05-27

# Checking original sequences with FastQC

We used: <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>)

# Use seqtk to clean reads (recommended by CBW)

## need to download to server:

``` bash
git clone git://github.com/lh3/seqtk.git
#./seqtk/seqtk
# put in same directory as samples
```

## make a sample list

needed to ls \*fastq, but narrow window so each sample appears on a line
copied pasted into vim sample\_list - need to press i to insert, press
ESC, - then colon (:) and x to save (or q\! to not save), and quit

## using seqtk on all samples

made vim trim\_samples.sh

Make a for loop:

``` bash
# make a shell script in vim (trim_samples.sh)
i
--
while read name;
do ./seqtk/seqtk trimfq -q 0.05 "$name" > trim_"$name";

done < sample_list
--
ESC 
:x

# saved, and then ran
nohup bash trim_samples.sh &
  
# make a list of trimmed samples using vim (just like above)
trim_sample_list  
```

# Reference sequence prep

## sequence download, then unzip

``` bash
# KH35c
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/197/105/GCA_002197105.1_ASM219710v1/GCA_002197105.1_ASM219710v1_genomic.fna.gz

# unzip
gunzip GCA_002197105.1_ASM219710v1_genomic.fna.gz
```

## index reference using bwa

``` bash
/gran1/apps/bwa-0.7.13/bwa index GCA_002197105.1_ASM219710v1_genomic.fna
# creates a bunch of .fna files w/different extensions
```

## extract reference w/ samtools

``` bash
/gran1/apps/samtools-1.3/samtools faidx GCA_002197105.1_ASM219710v1_genomic.fna
# creates .fai file
```

## make the .dict file using picard

``` bash
java -jar /gran1/rebecca.batstone/Alignment_to_VCF/picard/build/libs/picard.jar CreateSequenceDictionary REFERENCE=GCA_002197105.1_ASM219710v1_genomic.fna OUTPUT=GCA_002197105.1_ASM219710v1_genomic.fna.dict
# makes a .dict file
```

# Alignment: GATK best practices workflow

I need to transfer the original reads from capsicum to grandiflora. To
zip, I used: nohup zip -r FASTQ\_sequences.zip FASTQ\_WGS\_2017 \>
zipped\_seq.out &

Ended up just wget’ing the original reads, re-doing seqtk. Now everthing
is in align\_KH35c\_08Mar2018

We begin by mapping the sequence reads to the reference genome to
produce a file in SAM/BAM format sorted by coordinate. Next, we mark
duplicates to mitigate biases introduced by data generation steps such
as PCR amplification. Finally, we recalibrate the base quality scores,
because the variant calling algorithms rely heavily on the quality
scores assigned to the individual base calls in each sequence read.

## Align reads to reference

map reads to reference (mem works better if sequence quality is good)

“The variant discovery workflow depends on having sequence data in the
form of reads that are aligned to a reference genome. So the very first
step is of course to map your reads to the reference to produce a file
in SAM/BAM format.”

Note about modifying read names e.g.,
trim\_A01-1021-11\_171006\_NextSeq\_R1.fastq
n=\({name%_Next*} # n becomes trim_A01-1021-11_171006 n=\){name%\_R\*}
\# n becomes trim\_A01-1021-11\_171006\_NextSeq

``` bash
# make a shell script using vim (bwa_mem.sh)
i
--  
while read name;
do /gran1/apps/bwa-0.7.13/bwa mem -t 4 -M -R "@RG\tID:"${name%_171006*}"\tSM:"${name%_171006*}"" GCA_002197105.1_ASM219710v1_genomic.fna "${name%_R*}"_R1.fastq "${name%_R*}"_R2.fastq > "${name%_171006*}".sam;

done < trimmed_samples.list
--
ESC
:x

# run (in background)
nohup bash bwa_mem.sh > bwa_mem.out &
  
# run bwa for WSM1022
nohup /gran1/apps/bwa-0.7.13/bwa mem -t 4 -M -R "@RG\tID:WSM1022\tSM:WSM1022" GCA_002197105.1_ASM219710v1_genomic.fna GCA_000510665.1_ASM51066v1_genomic.fna > WSM1022.sam &
  
# make a list of samples using vim
sam_samples.list  
```

## Convert sam to bam

``` bash
# make a shell script using vim (sam_to_bam.sh)
i
--
while read name;
do /gran1/apps/samtools-1.3/samtools view -huS -o "${name%.sam*}".bam "$name";

done < sam_samples.list # made list same way as previous
--
ESC
:x

# run (in background)
nohup bash sam_to_bam.sh > sam_to_bam.out &
  
# make a list of samples in vim
bam_samples.list
```

## Reorder using PICARD

``` bash
# make a shell script using vim (picard_reorder.sh)
i
--
while read name;
do java -jar /gran1/rebecca.batstone/Alignment_to_VCF/picard/build/libs/picard.jar ReorderSam R=GCA_002197105.1_ASM219710v1_genomic.fna I="$name" O="$name".reorder;

done < bam_samples.list # made list same way as previous
--
ESC
:x

# run (in background)
nohup bash picard_reorder.sh > picard_reorder.out &
  
# make a sample list using vim
reorder_samples.list
```

## Add or replace read groups using Picard

“Read group information is typically added during this step, but can
also be added or modified after mapping using Picard
AddOrReplaceReadGroups.”

Note, I forgot to do this here but did it later (just before
HaplotypeCaller) and ran into many errors. So I had to restart from this
point…

Also, you can do both sorted (SO=coordinate) and indexing
(creat\_index=TRUE) in one go, but I decided against it so I could make
sure things were working more frequently. Also, markdups should be done
before indexing anyways.

Notes about previous runs:

parts I took out: CREATE\_INDEX= True \# changed to false, index after
(see below) SO= coordinate \# do after in another step

Error: for sample F06-1022-13 in indexing step, apparently the aligned
reads were too large. But I didn’t get this error before, so I thought
to just index after. The actual error is below:

Exception when processing alignment for BAM index 2/2 151b aligned read

``` bash
# make a shell script (add_read_groups.sh), and then bash it
vim add_read_groups.sh
i
--
while read name;
do java -jar /gran1/rebecca.batstone/Alignment_to_VCF/picard/build/libs/picard.jar AddOrReplaceReadGroups I="$name" O="$name".RG RGLB=Rhizo_evo RGPU=1 RGPL=Illumina RGSM="$name" 

done < reorder_samples.list
--
ESC
:x

# run shell script in background
nohup bash add_RG.sh > add_RG.out &
  
# after it runs, make a list of samples w/ read groups added
vim RG_samples.list
PASTE all names
ESC
:x
```

## Sorting using Picard

“Once you have mapped the reads, you’ll need to make sure they are
sorted in the proper order (by coordinate).”

``` bash
# make a shell script using vim (sort_sam.sh)
i
--
while read name;
do java -jar /gran1/rebecca.batstone/Alignment_to_VCF/picard/build/libs/picard.jar SortSam I="$name" O="$name".sorted SO=coordinate;

done < RG_samples.list 
--
ESC
:x

# run (in background)
nohup bash sort_sam.sh > sort_sam.sh &
  
# make a list of samples using vim
sorted_samples.list
```

## Mark duplicates

“Once your data has been mapped to the reference genome, you can proceed
to mark duplicates. The idea here is that during the sequencing process,
the same DNA fragments may be sequenced several times. The resulting
duplicate reads are not informative and should not be counted as
additional evidence for or against a putative variant. The duplicate
marking process (sometimes called *dedupping* in bioinformatics slang)
does not remove the reads, but identifies them as duplicates by adding a
flag in the read’s SAM record. Most GATK tools will then ignore these
duplicate reads by default, through the internal application of a read
filter.”

\[Marking reads that look suspiciously similar\]

another error (after running w/ create\_index=false): Write error;
BinaryCodec in writemode; streamed file (filename not available) This
error was b/c capsicum was full\!

So, I started using grandiflora downstream of the sort step. The code
remains the same except for specifying the program paths (see above)

``` bash
# make a shell script (markdups.sh)
i
--
while read name;
do java -jar /gran1/rebecca.batstone/Alignment_to_VCF/picard/build/libs/picard.jar MarkDuplicates I="$name" O="$name".markdups M=duplicates;

done < sorted_samples.list 
--

# run (in background)
nohup bash markdups.sh > markdups.out &
  
# make a list of samples using vim (note the period rather than underscore)
markdups_samples.list
```

## Index the bam (after marked dups)

``` bash
# make a shell script using vim (index.sh)
i
--
while read name;
do java -jar /gran1/rebecca.batstone/Alignment_to_VCF/picard/build/libs/picard.jar BuildBamIndex I="$name";

done < markdups_samples.list
--
  
# run (in background)
nohup bash index.sh > index.out &
```

## Base recalibration

NOTE: I didn’t complete this step “Variant calling algorithms rely
heavily on the quality scores assigned to the individual base calls in
each sequence read. These scores are per-base estimates of error emitted
by the sequencing machines. Unfortunately the scores produced by the
machines are subject to various sources of systematic technical error,
leading to over- or under-estimated base quality scores in the data.
Base quality score recalibration (BQSR) is a process in which we apply
machine learning to model these errors empirically and adjust the
quality scores accordingly. This allows us to get more accurate base
qualities, which in turn improves the accuracy of our variant calls.”

# GATK variant discovery

This workflow involves running HaplotypeCaller on each sample separately
in GVCF mode, to produce an intermediate file format called GVCF (for
Genomic VCF), which is described in detail in the documentation. The
GVCFs of multiple samples are then run through a joint genotyping step
to produce a multi-sample VCF callset, which can then be filtered to
balance sensitivity and specificity as desired. In the final analysis,
this workflow produces results equivalent to traditional joint calling,
in which all samples are given simultaneously to the variant caller, but
it scales much better and resolves the so-called N+1 problem. Note that
this is equally applicable to small cohorts or even single samples.

## Haplotype caller

for many individuals, use this. this makes a file calling all the SNPs
relative to the reference

Notes: \* first step, on individual bams, then genotypeGVCF on all\! \*
takes a really long time per sample (30 hours +). \* Fastest way to
process is to divide samples up into equal-sized groups (in my case, 5-6
samples each) and run up to ten haplotype callers at a time. \* For me,
took about one week to complete. \* needs .bam extention, or haplotype
caller won’t run \* previous calls were for diploid organisms, used new
gatk HaplotypeCaller to call SNPs with haploid settings

``` bash
# example renaming:
Rename all *.markdups to *.markdups.bam
for f in *.markdups; do 
mv -- "$f" "${f%.markdups}.markdups.bam"
done

# e.g., hap_call_1-5.sh:
while read name;
do $GATK HaplotypeCaller -R GCA_002197105.1_ASM219710v1_genomic.fna -I "$name" --dont-use-soft-clipped-bases TRUE -ploidy 1 -O "$name".g.vcf -ERC GVCF;

done < samples_1-5.list

nohup bash hap_call_1-5.sh > hap_call_1-5.out &
  
# ran all 10 .sh files at a time. Took ~ 3 days
```

## Run Combine GVCFs

``` bash
nohup $GATK CombineGVCFs -R GCA_002197105.1_ASM219710v1_genomic.fna -V gVCFs.list -O comb_KH35c.vcf > combine_gVCFs.out &
```

## GenotypeGVCFs (on grandiflora)

“At this step, we gather all the per-sample GVCFs (or combined GVCFs if
we are working with large numbers of samples) and pass them all together
to the joint genotyping tool, GenotypeGVCFs. This produces a set of
joint-called SNP and indel calls ready for filtering. This cohort-wide
analysis empowers sensitive detection of variants even at difficult
sites, and produces a squared-off matrix of genotypes that provides
information about all sites of interest in all samples considered, which
is important for many downstream analyses.”

This step runs very fast and can be rerun at any point when samples are
added to the cohort, thereby solving the so-called N+1 problem.

\[this makes a file containing a concatination of the .vcf’s made in the
previous step, so that you can compare between
samples.\]

``` bash
nohup $GATK GenotypeGVCFs -R GCA_002197105.1_ASM219710v1_genomic.fna -V comb_KH35c.vcf -ploidy 1 -O KH35c.vcf -stand-call-conf 30 > genotype_gVCFs.out &
```

## Rename samples in vcf (to make it easier downstream)

``` bash
bcftools reheader -s, --samples VCFs.list KH35c.vcf -o KH35c_rn.vcf
bcftools query -l KH35c_rn.vcf
```

## Variant quality score recalibration

NOTE: didn’t complete this step “The best way to filter the raw variant
callset is to use variant quality score recalibration (VQSR), which uses
machine learning to identify annotation profiles of variants that are
likely to be real, and assigns a VQSLOD score to each variant that is
much more reliable than the QUAL score calculated by the caller. In the
first step of this two-step process, the program builds a model based on
training variants, then applies that model to the data to assign a
well-calibrated probability to each variant call. We can then use this
variant quality score in the second step to filter the raw call set,
thus producing a subset of calls with our desired level of quality,
fine-tuned to balance specificity and sensitivity.”

“The downside of how variant recalibration works is that the algorithm
requires high-quality sets of known variants to use as training and
truth resources, which for many organisms are not yet available. It also
requires quite a lot of data in order to learn the profiles of good
vs. bad variants, so it can be difficult or even impossible to use on
small datasets that involve only one or a few samples, on targeted
sequencing data, on RNAseq, and on non-model organisms. If for any of
these reasons you find that you cannot perform variant recalibration on
your data (after having tried the workarounds that we recommend, where
applicable), you will need to use hard-filtering instead. This consists
of setting flat thresholds for specific annotations and applying them to
all variants equally.”

Used the “hard-filtering option” instead (see Variant\_filtering.md)