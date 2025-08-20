import pandas as pd

# --- Configuracion del Pipeline ---
configfile: "config/config.yaml"
SAMPLES = pd.read_csv(config['samples'], sep='\t').set_index('sample', drop=False)

REF_GENOME = "/opt/db/genome.fasta"
SNPEFF_DB = "Lbraziliensis_2019_manual"

# --- Regla Final ---
rule all:
    input:
        expand("results/annotated/{sample}.ann.vcf", sample=SAMPLES.index),
        "results/qc/multiqc_report.html"

# --- QC y Trimming ---
rule fastqc:
    input:
        r1="data/raw_fastq/{sample}_R1.fastq.gz",
        r2="data/raw_fastq/{sample}_R2.fastq.gz"
    output:
        zip_r1="results/qc/{sample}_R1_fastqc.zip",
        zip_r2="results/qc/{sample}_R2_fastqc.zip"
    params:
        outdir="results/qc"
    log: "results/logs/fastqc/{sample}.log"
    shell: "fastqc {input.r1} {input.r2} -o {params.outdir} --quiet &> {log}"

rule multiqc:
    input:
        expand("results/qc/{sample}_{read}_fastqc.zip", sample=SAMPLES.index, read=["R1", "R2"])
    output:
        "results/qc/multiqc_report.html"
    log: "results/logs/multiqc.log"
    shell: "multiqc results/qc -o results/qc -n multiqc_report.html --quiet &> {log}"

rule fastp:
    input:
        r1="data/raw_fastq/{sample}_R1.fastq.gz",
        r2="data/raw_fastq/{sample}_R2.fastq.gz"
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz"
    log: "results/logs/fastp/{sample}.log"
    threads: 8
    shell: "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} --thread {threads} -h /dev/null -j /dev/null &> {log}"

# --- Indexado de Referencia ---
rule bwa_index:
    input: REF_GENOME
    output: touch(REF_GENOME + ".bwa_indexed")
    log: "results/logs/bwa_index.log"
    shell: "bwa index {input} &> {log}"

rule samtools_faidx:
    input: REF_GENOME
    output: REF_GENOME + ".fai"
    log: "results/logs/samtools_faidx.log"
    shell: "samtools faidx {input} &> {log}"

rule gatk_create_dictionary:
    input: REF_GENOME
    output: REF_GENOME.replace(".fasta", ".dict")
    log: "results/logs/gatk_create_dictionary.log"
    resources:
        mem_gb=4
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{resources.mem_gb}g"
    shell: "gatk --java-options '{params.java_opts}' CreateSequenceDictionary -R {input} -O {output} &> {log}"

# --- Alineamiento y Post-procesamiento ---
rule bwa_mem_sort:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        ref=REF_GENOME,
        indexed=REF_GENOME + ".bwa_indexed"
    output:
        bam="results/aligned/{sample}.sorted.bam"
    params:
        read_group=r"'@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'"
    log: "results/logs/bwa_mem/{sample}.log"
    threads: 8
    shell: "bwa mem -t {threads} -R {params.read_group} {input.ref} {input.r1} {input.r2} | samtools view -bS - | samtools sort -o {output.bam} &> {log}"

rule mark_duplicates:
    input: "results/aligned/{sample}.sorted.bam"
    output: "results/aligned/{sample}.dedup.bam"
    log: "results/logs/mark_duplicates/{sample}.log"
    resources:
        mem_gb=8
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{resources.mem_gb}g"
    shell: "gatk --java-options '{params.java_opts}' MarkDuplicates -I {input} -O {output} -M /dev/null &> {log}"

rule samtools_index:
    input: "results/aligned/{sample}.dedup.bam"
    output: "results/aligned/{sample}.dedup.bam.bai"
    log: "results/logs/samtools_index/{sample}.log"
    shell: "samtools index {input} &> {log}"

# --- Llamada de Variantes ---
rule haplotype_caller:
    input:
        bam="results/aligned/{sample}.dedup.bam",
        bai="results/aligned/{sample}.dedup.bam.bai",
        ref=REF_GENOME,
        fai=REF_GENOME + ".fai",
        dict=REF_GENOME.replace(".fasta", ".dict")
    output:
        vcf="results/variants/{sample}.vcf.gz"
    resources:
        mem_gb=config.get('max_mem_gb', 16)
    params:
        java_opts=lambda wildcards, resources: f"-Xmx{resources.mem_gb}g"
    log: "results/logs/haplotype_caller/{sample}.log"
    shell: "gatk --java-options '{params.java_opts}' HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} &> {log}"

# --- Anotacion ---
rule snpeff_annotate:
    input:
        vcf="results/variants/{sample}.vcf.gz"
    output:
        vcf="results/annotated/{sample}.ann.vcf",
        html="results/annotated/{sample}.snpeff.html"
    params:
        db=SNPEFF_DB
    log: "results/logs/snpeff/{sample}.log"
    shell: "snpEff ann -v {params.db} {input.vcf} > {output.vcf} -stats {output.html} 2> {log}"
