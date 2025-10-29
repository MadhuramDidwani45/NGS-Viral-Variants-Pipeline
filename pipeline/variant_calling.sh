#!/bin/bash

set -e  

# Configuration
REFERENCE="reference.fa"
THREADS=12

# Sample list
SAMPLES=(SRR33784444 SRR33784445 SRR33784446 SRR33784447 SRR33784448 SRR33784449)

echo "Starting HIV-1 variant calling pipeline..."
date

# Step 0: Create required directories
echo "Step 0: Creating required directories..."
mkdir -p trimmed_reads
mkdir -p bam_files
mkdir -p vcf_files
mkdir -p reports
mkdir -p raw_fastq
mkdir -p Fastqc_results


# Step 0.1: Download raw FASTQ files
echo "Step 0.1: Downloading raw FASTQ files..."
for SAMPLE in "${SAMPLES[@]}"; do
    echo "  Downloading ${SAMPLE}..."
    fastq-dump --split-files --gzip --outdir raw_fastq ${SAMPLE}
done

# Step 0.2: Quality Control FASTQC
echo "Step 0.2: Quality COntrol with Fastqc..."
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Fastqc ${SAMPLE}..."
    fastqc ${SAMPLE} -o Fastqc_results/
done

# Step 0.3: Quality trimming with Trimmomatic
echo "Step 0.3: Trimming reads with Trimmomatic..."
for SAMPLE in "${SAMPLES[@]}"; do
    echo "  Trimming ${SAMPLE}..."
    trimmomatic PE -threads ${THREADS} \
        raw_fastq/${SAMPLE}_1.fastq.gz \
        raw_fastq/${SAMPLE}_2.fastq.gz \
        trimmed_reads/${SAMPLE}_1_paired.fastq.gz \
        trimmed_reads/${SAMPLE}_1_unpaired.fastq.gz \
        trimmed_reads/${SAMPLE}_2_paired.fastq.gz \
        trimmed_reads/${SAMPLE}_2_unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# Step 1: Download reference genome
echo "Step 1: Downloading HIV-1 reference genome..."
wget --quiet --show-progress \
    --output-document=${REFERENCE} \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001802.1&rettype=fasta"

# Step 2: Index reference genome
echo "Step 2: Indexing reference genome..."
bwa index ${REFERENCE}
samtools faidx ${REFERENCE}
gatk CreateSequenceDictionary -R ${REFERENCE}

# Step 3: Alignment with BWA-MEM
echo "Step 3: Aligning reads with BWA-MEM..."
for SAMPLE in "${SAMPLES[@]}"; do
    echo "  Processing ${SAMPLE}..."
    bwa mem -t ${THREADS} -M \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
        ${REFERENCE} \
        trimmed_reads/${SAMPLE}_1_paired.fastq.gz \
        trimmed_reads/${SAMPLE}_2_paired.fastq.gz > bam_files/${SAMPLE}_aligned.sam
done

# Step 4: Convert SAM to BAM and sort
echo "Step 4: Converting to BAM and sorting..."
for SAMPLE in "${SAMPLES[@]}"; do
    samtools view -bS bam_files/${SAMPLE}_aligned.sam | \
        samtools sort -o bam_files/${SAMPLE}_sorted.bam
    rm bam_files/${SAMPLE}_aligned.sam
done

# Step 5: Index BAM files
echo "Step 5: Indexing BAM files..."
for SAMPLE in "${SAMPLES[@]}"; do
    samtools index bam_files/${SAMPLE}_sorted.bam
done

# Step 6: Generate alignment statistics
echo "Step 6: Generating alignment statistics..."
for SAMPLE in "${SAMPLES[@]}"; do
    samtools flagstat bam_files/${SAMPLE}_sorted.bam > reports/${SAMPLE}_alignment_stats.txt
done

# Step 7: Calculate coverage
echo "Step 7: Calculating coverage..."
for SAMPLE in "${SAMPLES[@]}"; do
    samtools depth bam_files/${SAMPLE}_sorted.bam > reports/${SAMPLE}_coverage.txt
    awk '{sum+=$3; count++} END {print "Average coverage:", sum/count}' \
        reports/${SAMPLE}_coverage.txt > reports/${SAMPLE}_coverage_summary.txt
done

# Step 8: Mark duplicates
echo "Step 8: Marking duplicates..."
for SAMPLE in "${SAMPLES[@]}"; do
    gatk MarkDuplicates \
        -I bam_files/${SAMPLE}_sorted.bam \
        -O bam_files/${SAMPLE}_dedup.bam \
        -M reports/${SAMPLE}_dup_metrics.txt \
        --REMOVE_DUPLICATES false \
        --CREATE_INDEX true
done

# Step 9: Variant calling with Mutect2
echo "Step 9: Calling variants with Mutect2..."
for SAMPLE in "${SAMPLES[@]}"; do
    echo "  Calling variants for ${SAMPLE}..."
    gatk Mutect2 \
        -R ${REFERENCE} \
        -I bam_files/${SAMPLE}_dedup.bam \
        -O vcf_files/${SAMPLE}_raw.vcf \
        --native-pair-hmm-threads ${THREADS} \
        --max-reads-per-alignment-start 0 \
        --min-base-quality-score 20 \
        --callable-depth 10
done

# Step 10: Filter variants
echo "Step 10: Filtering variants..."
for SAMPLE in "${SAMPLES[@]}"; do
    gatk FilterMutectCalls \
        -R ${REFERENCE} \
        -V vcf_files/${SAMPLE}_raw.vcf \
        -O vcf_files/${SAMPLE}_filtered.vcf
done

# Step 11: Select PASS variants only
echo "Step 11: Selecting PASS variants..."
for SAMPLE in "${SAMPLES[@]}"; do
    gatk SelectVariants \
        -R ${REFERENCE} \
        -V vcf_files/${SAMPLE}_filtered.vcf \
        -O vcf_files/${SAMPLE}_pass.vcf \
        --exclude-filtered true
done

# Step 12: Generate variant statistics
echo "Step 12: Generating variant statistics..."
for SAMPLE in "${SAMPLES[@]}"; do
    bcftools stats vcf_files/${SAMPLE}_pass.vcf > reports/${SAMPLE}_variant_stats.txt
done

# Step 13: Merge all samples
echo "Step 13: Merging all samples..."
for f in vcf_files/*_pass.vcf; do
    bgzip -c "$f" > "${f}.gz"
done

for f in vcf_files/*_pass.vcf.gz; do
    tabix -p vcf "$f"
done

bcftools merge vcf_files/*_pass.vcf.gz -o vcf_files/merged_all_samples.vcf

echo "Pipeline completed successfully!"
date