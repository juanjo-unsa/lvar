# ==============================================================================
# Snakefile para LVAR: Analisis de Variantes en Leishmania
# ==============================================================================
#
# Autor: Juanjo (con ayuda de IA)
# Descripcion: Pipeline para el analisis de WGS desde lecturas crudas
#              hasta la anotacion de variantes, usando Snakemake y Docker.
#
# ==============================================================================

import pandas as pd

# --- Carga de Configuracion ---
configfile: "config/config.yaml"
SAMPLES = pd.read_csv(config['samples'], sep='\t').set_index('sample', drop=False)

# --- Rutas y Parametros Globales (dentro del contenedor) ---
# Usamos la base de datos que construimos manualmente en el Dockerfile
SNPEFF_DB = "Lbraziliensis_2019_manual"
# Ruta al genoma que preparamos en el Dockerfile
REF_GENOME = "/opt/db/genome.fasta"

# ==============================================================================
# REGLA FINAL (TARGET)
# ==============================================================================
rule all:
    input:
        expand("results/annotated/{sample}.ann.vcf", sample=SAMPLES.index),
        "results/qc/multiqc_report.html"

# ==============================================================================
# PASO 1: QC
# ==============================================================================
rule fastqc:
    input: r1="data/raw_fastq/{sample}_R1.fastq.gz", r2="data/raw_fastq/{sample}_R2.fastq.gz"
    output: html_r1="results/qc/{sample}_R1_fastqc.html", zip_r1="results/qc/{sample}_R1_fastqc.zip", html_r2="results/qc/{sample}_R2_fastqc.html", zip_r2="results/qc/{sample}_R2_fastqc.zip"
    params: outdir="results/qc"
    log: "results/logs/fastqc/{sample}.log"
    shell: "fastqc {input.r1} {input.r2} -o {params.outdir} > {log} 2>&1"

rule multiqc:
    input: expand("results/qc/{sample}_{read}_fastqc.zip", sample=SAMPLES.index, read=["R1", "R2"])
    output: "results/qc/multiqc_report.html"
    log: "results/logs/multiqc.log"
    shell: "multiqc results/qc -o results/qc -n multiqc_report.html > {log} 2>&1"

# ==============================================================================
# PASO 2: TRIMMING
# ==============================================================================
rule fastp:
    input: r1="data/raw_fastq/{sample}_R1.fastq.gz", r2="data/raw_fastq/{sample}_R2.fastq.gz"
    output: r1="results/trimmed/{sample}_R1.fastq.gz", r2="results/trimmed/{sample}_R2.fastq.gz", html="results/qc/{sample}.fastp.html"
    log: "results/logs/fastp/{sample}.log"
    threads: 8
    shell: "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -h {output.html} --thread {threads} > {log} 2>&1"

# ==============================================================================
# PASO 3: ALINEAMIENTO
# ==============================================================================
rule bwa_index:
    input: REF_GENOME
    output: touch(REF_GENOME + ".bwa_indexed")
    log: "results/logs/bwa_index.log"
    shell: "bwa index {input} > {log} 2>&1"

rule bwa_mem_sort:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        ref=REF_GENOME,
        indexed=REF_GENOME + ".bwa_indexed"
    output: bam="results/aligned/{sample}.sorted.bam"
    params: read_group=r"'@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'"
    log: "results/logs/bwa_mem/{sample}.log"
    threads: 8
    shell: "bwa mem -t {threads} -R {params.read_group} {input.ref} {input.r1} {input.r2} | samtools view -bS - | samtools sort -o {output.bam} > {log} 2>&1"

# ==============================================================================
# PASO 4: POST-PROCESAMIENTO
# ==============================================================================
rule mark_duplicates:
    input: "results/aligned/{sample}.sorted.bam"
    output: bam="results/aligned/{sample}.dedup.bam", metrics="results/aligned/{sample}.dedup.metrics"
    log: "results/logs/mark_duplicates/{sample}.log"
    shell: "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} > {log} 2>&1"

rule samtools_index:
    input: "results/aligned/{sample}.dedup.bam"
    output: "results/aligned/{sample}.dedup.bam.bai"
    log: "results/logs/samtools_index/{sample}.log"
    shell: "samtools index {input} > {log} 2>&1"

# ==============================================================================
# PASO 5: LLAMADA DE VARIANTES
# ==============================================================================
rule haplotype_caller:
    input:
        bam="results/aligned/{sample}.dedup.bam",
        bai="results/aligned/{sample}.dedup.bam.bai",
        ref=REF_GENOME
    output:
        vcf="results/variants/{sample}.vcf.gz"
    params:
        java_opts=f"-Xmx{config.get('gatk_ram_gb', 4)}g"
    log: "results/logs/haplotype_caller/{sample}.log"
    shell: "gatk --java-options '{params.java_opts}' HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} > {log} 2>&1"

# ==============================================================================
# PASO 6: ANOTACION CON SnpEff
# ==============================================================================
rule snpeff_annotate:
    input:
        vcf="results/variants/{sample}.vcf.gz"
    output:
        vcf="results/annotated/{sample}.ann.vcf",
        html="results/annotated/{sample}.snpeff.html"
    params:
        db=SNPEFF_DB
    log: "results/logs/snpeff/{sample}.log"
    shell:
        """
        snpEff ann -v {params.db} {input.vcf} > {output.vcf} \
        -stats {output.html} 2> {log}
        """
