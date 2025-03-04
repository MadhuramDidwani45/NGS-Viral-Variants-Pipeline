#!/bin/bash
set -e

LOGFILE="variant_calling.log"
exec > >(tee -a "$LOGFILE") 2>&1

# Define a logging function with a timestamp.
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log "Starting variant calling pipeline..."

# Prefetch and fastq-dump steps.
log "Prefetching SRR32313970..."
prefetch SRR32313970 --progress
log "Prefetch complete."

log "Running fastq-dump to split FASTQ files..."
fastq-dump --split-files SRR32313970/
log "fastq-dump complete."

# Run FastQC on raw FASTQ files.
log "Running FastQC on raw FASTQ files..."
fastqc SRR32313970_1.fastq SRR32313970_2.fastq
log "FastQC on raw FASTQ files complete."

# Create directory for trimmed reads and run Trimmomatic.
log "Creating directory for trimmed reads and running Trimmomatic..."
mkdir  Trimmed
cd Trimmed
trimmomatic PE ../SRR32313970_1.fastq ../SRR32313970_2.fastq \
    SRR32313970_1_paired.fastq SRR32313970_1_unpaired.fastq \
    SRR32313970_2_paired.fastq SRR32313970_2_unpaired.fastq \
    SLIDINGWINDOW:4:20 MINLEN:50
log "Trimmomatic processing complete."

# Run FastQC on trimmed files.
log "Running FastQC on trimmed FASTQ files..."
fastqc SRR32313970_1_paired.fastq SRR32313970_2_paired.fastq SRR32313970_1_unpaired.fastq SRR32313970_2_unpaired.fastq
log "FastQC on trimmed FASTQ files complete."
cd ..

# Alignment with Bowtie2.
log "Changing to Trimmed directory for alignment..."
cd Trimmed
log "Build na index for the reference genome"
bowtie2-build GRCh38.primary_assembly.genome.fa index
log "Running Bowtie2 alignment..."
bowtie2 --no-unal -p 2 -x /home/madhuram9011/Downloads/variant_call_pipe/trimmed_reads/index \
    -1 SRR32313970_1_paired.fastq -2 SRR32313970_2_paired.fastq -S alignment.sam
log "Bowtie2 alignment complete."
mv alignment.sam ../
cd ..

# Convert SAM to BAM, sort, and index using SAMtools.
log "Converting SAM to BAM..."
samtools view -Sb -o alignment.bam alignment.sam
log "Sorting BAM file..."
samtools sort -O bam -o sorted.bam alignment.bam
log "Indexing sorted BAM file..."
samtools index sorted.bam

# Add or replace read groups with GATK.
log "Adding or replacing read groups with GATK..."
gatk AddOrReplaceReadGroups \
    -I sorted.bam \
    -O sorted_rg.bam \
    --RGID 1 \
    --RGLB lib1 \
    --RGPL ILLUMINA \
    --RGPU unit1 \
    --RGSM SampleName
log "Read groups added/replaced."

# Mark duplicates using GATK MarkDuplicates.
log "Marking duplicates with GATK..."
gatk MarkDuplicates \
    -I sorted_rg.bam \
    -R GRCh38.primary_assembly.genome.fa  \
    -M metrics.txt \
    -O unique_reads.bam
log "Duplicates marked."

# Rescore the bases using GATK BaseRecalibrator.
log "Performing BaseRecalibration with GATK BaseRecalibrator..."
gatk BaseRecalibrator \
    -R GRCh38.primary_assembly.genome.fa  \
    -I unique_reads.bam \
    --known-sites Mills_and_1000G_gold_standard.indels.hg38.renamed.vcf.gz \
    -O recal_data.table
log "BaseRecalibration complete."

log "Applying BQSR with GATK ApplyBQSR..."
gatk ApplyBQSR \
    -R GRCh38.primary_assembly.genome.fa  \
    -I unique_reads.bam \
    --bqsr-recal-file recal_data.table \
    -O recalibrated.bam
log "BQSR applied."

# Call variants using GATK HaplotypeCaller.
log "Calling variants with GATK HaplotypeCaller..."
gatk HaplotypeCaller \
    -R GRCh38.primary_assembly.genome.fa  \
    -I recalibrated.bam \
    -O output.vcf
log "Variant calling complete."

# Filter variants using GATK VariantFiltration.
log "Filtering variants with GATK VariantFiltration..."
gatk VariantFiltration \
    -R GRCh38.primary_assembly.genome.fa  \
    -V output.vcf \
    -O filtered_output.vcf \
    --filter-expression "QD < 2.0" \
    --filter-name "LowQD" \
    --filter-expression "FS > 60.0" \
    --filter-name "HighFS"
log "Variant filtration complete."

# Annotate variants using GATK VariantAnnotator.
log "Annotating variants with GATK VariantAnnotator..."
gatk VariantAnnotator \
    -R GRCh38.primary_assembly.genome.fa \
    -V filtered_output.vcf \
    -I recalibrated.bam \
    -O output.annotated.vcf \
    -A Coverage \
    -A QualByDepth \
    -A MappingQualityRankSumTest
log "Variant annotation complete."


log "Variant calling pipeline finished successfully."
