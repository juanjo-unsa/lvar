# Snakefile para análisis de variantes en Leishmania braziliensis
# ---------------------------------------------------------------

import pandas as pd

# Cargar la configuración y la lista de muestras
configfile: "config/config.yaml"
SAMPLES = pd.read_csv(config['samples'], sep='\t').set_index('sample', drop=False)

# --- REGLA FINAL: DEFINE LOS ARCHIVOS FINALES DESEADOS ---
# Snakemake intentará generar estos archivos, ejecutando las reglas necesarias.
rule all:
    input:
        # Anotación final para cada muestra
        expand("results/annotated/{sample}.ann.vcf", sample=SAMPLES.index),
        # Reporte de calidad agregado
        "results/qc/multiqc_report.html"

# --- PASO 1: Control de calidad con FastQC y MultiQC ---
rule fastqc:
    input:
        r1="data/raw_fastq/{sample}_R1.fastq.gz",
        r2="data/raw_fastq/{sample}_R2.fastq.gz"
    output:
        html_r1="results/qc/{sample}_R1_fastqc.html",
        zip_r1="results/qc/{sample}_R1_fastqc.zip",
        html_r2="results/qc/{sample}_R2_fastqc.html",
        zip_r2="results/qc/{sample}_R2_fastqc.zip"
    params:
        outdir="results/qc"
    log:
        "results/logs/fastqc/{sample}.log"
    container:
        "docker://staphb/fastqc:latest"
    shell:
        "fastqc {input.r1} {input.r2} -o {params.outdir} > {log} 2>&1"

rule multiqc:
    input:
        expand("results/qc/{sample}_{read}_fastqc.zip", sample=SAMPLES.index, read=["R1", "R2"])
    output:
        "results/qc/multiqc_report.html"
    log:
        "results/logs/multiqc.log"
    container:
        "docker://ewels/multiqc:latest"
    shell:
        "multiqc results/qc -o results/qc -n multiqc_report.html > {log} 2>&1"

# --- PASO 2: Limpieza de lecturas con fastp ---
rule fastp:
    input:
        r1="data/raw_fastq/{sample}_R1.fastq.gz",
        r2="data/raw_fastq/{sample}_R2.fastq.gz"
    output:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        html="results/qc/{sample}.fastp.html"
    log:
        "results/logs/fastp/{sample}.log"
    threads: 8
    container:
        "docker://quay.io/biocontainers/fastp:0.23.2--h79da9fb_0"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-h {output.html} --thread {threads} > {log} 2>&1"

# --- PASO 3: Indexar el genoma de referencia para BWA ---
rule bwa_index:
    input:
        config["ref_genome"]
    output:
        touch(config["ref_genome"] + ".bwa_indexed")
    log:
        "results/logs/bwa_index.log"
    container:
        "docker://biocontainers/bwa:v0.7.17_cv1"
    shell:
        "bwa index {input} > {log} 2>&1"

# --- PASO 4: Alineamiento con BWA-MEM y ordenamiento con Samtools ---
rule bwa_mem_sort:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        ref=config["ref_genome"],
        indexed=config["ref_genome"] + ".bwa_indexed"
    output:
        bam="results/aligned/{sample}.sorted.bam"
    params:
        read_group=r"'@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'"
    log:
        "results/logs/bwa_mem/{sample}.log"
    threads: 8
    container:
        "docker://biocontainers/bwa:v0.7.17_cv1"
    shell:
        "bwa mem -t {threads} -R {params.read_group} {input.ref} {input.r1} {input.r2} | "
        "samtools view -bS - | "
        "samtools sort -o {output.bam} > {log} 2>&1"

# --- PASO 5: Marcar duplicados con GATK ---
rule mark_duplicates:
    input:
        "results/aligned/{sample}.sorted.bam"
    output:
        bam="results/aligned/{sample}.dedup.bam",
        metrics="results/aligned/{sample}.dedup.metrics"
    log:
        "results/logs/mark_duplicates/{sample}.log"
    container:
        "docker://broadinstitute/gatk:4.2.6.1"
    shell:
        "gatk MarkDuplicates -I {input} -O {output.bam} -M {output.metrics} > {log} 2>&1"

# --- PASO 6: Indexar el BAM final ---
rule samtools_index:
    input:
        "results/aligned/{sample}.dedup.bam"
    output:
        "results/aligned/{sample}.dedup.bam.bai"
    log:
        "results/logs/samtools_index/{sample}.log"
    container:
        "docker://biocontainers/samtools:v1.9-4-deb_cv1"
    shell:
        "samtools index {input} > {log} 2>&1"

# --- PASO 7: Llamada de variantes con GATK HaplotypeCaller ---
rule haplotype_caller:
    input:
        bam="results/aligned/{sample}.dedup.bam",
        bai="results/aligned/{sample}.dedup.bam.bai",
        ref=config["ref_genome"]
    output:
        vcf="results/variants/{sample}.vcf.gz"
    log:
        "results/logs/haplotype_caller/{sample}.log"
    container:
        "docker://broadinstitute/gatk:4.2.6.1"
    shell:
        "gatk HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} > {log} 2>&1"

# --- PASO 8: Anotación de variantes con SnpEff ---
# ¡ATENCIÓN! Este paso requiere que tengas una base de datos de SnpEff para L. braziliensis.
# Consulta la documentación de SnpEff para crearla a partir de un archivo GFF y el FASTA de referencia.
# Comando general: java -jar snpEff.jar build -gff3 -v <nombre_db>
rule snpeff_annotate:
    input:
        vcf="results/variants/{sample}.vcf.gz",
        ref=config["ref_genome"] # solo como dependencia
    output:
        vcf="results/annotated/{sample}.ann.vcf",
        html="results/annotated/{sample}.snpeff.html"
    params:
        db=config["snpeff"]["database"]
    log:
        "results/logs/snpeff/{sample}.log"
    container:
        "docker://snpeff/snpeff:5.1d-1"
    shell:
        "snpEff ann -v {params.db} {input.vcf} > {output.vcf} "
        "-stats {output.html} 2> {log}"
