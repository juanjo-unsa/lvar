# ==============================================================================
# Snakefile para LVAR: Análisis de Variantes en Leishmania
# ==============================================================================
#
# Autor: Juanjo (con ayuda de IA)
# Descripción: Pipeline para el análisis de WGS desde lecturas crudas
#              hasta la anotación de variantes, usando Snakemake y Docker.
#
# ==============================================================================

import pandas as pd

# --- Carga de Configuración ---
# Carga los parámetros del pipeline desde el archivo YAML generado.
configfile: "config/config.yaml"

# Carga la lista de muestras desde el archivo TSV generado.
# Esto hace que el pipeline sea agnóstico al número y nombre de las muestras.
SAMPLES = pd.read_csv(config['samples'], sep='\t').set_index('sample', drop=False)

# ==============================================================================
# REGLA FINAL (TARGET)
# ==============================================================================
# Esta regla define los archivos finales que queremos que Snakemake genere.
# Snakemake trabajará hacia atrás desde aquí para ejecutar las reglas necesarias.
rule all:
    input:
        # Los VCF anotados para cada muestra.
        expand("results/annotated/{sample}.ann.vcf", sample=SAMPLES.index),
        # Un reporte de calidad agregado de todas las muestras.
        "results/qc/multiqc_report.html"

# ==============================================================================
# PASO 1: CONTROL DE CALIDAD (QC)
# ==============================================================================
# Ejecuta FastQC sobre las lecturas crudas y MultiQC para agregar los reportes.
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
        # Recolecta todos los archivos ZIP de FastQC para el reporte.
        expand("results/qc/{sample}_R1_fastqc.zip", sample=SAMPLES.index),
        expand("results/qc/{sample}_R2_fastqc.zip", sample=SAMPLES.index)
    output:
        "results/qc/multiqc_report.html"
    log:
        "results/logs/multiqc.log"
    container:
        "docker://ewels/multiqc:latest"
    shell:
        "multiqc results/qc -o results/qc -n multiqc_report.html > {log} 2>&1"

# ==============================================================================
# PASO 2: LIMPIEZA DE LECTURAS (TRIMMING)
# ==============================================================================
# Elimina adaptadores y bases de baja calidad con fastp.
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

# ==============================================================================
# PASO 3: ALINEAMIENTO
# ==============================================================================
# Primero, indexa el genoma de referencia para BWA (se ejecuta una sola vez).
rule bwa_index:
    input:
        config["ref_genome"]
    output:
        # Se crea un archivo "flag" para indicar que la indexación se completó.
        touch(config["ref_genome"] + ".bwa_indexed")
    log:
        "results/logs/bwa_index.log"
    container:
        "docker://biocontainers/bwa:v0.7.17_cv1"
    shell:
        "bwa index {input} > {log} 2>&1"

# Alinea las lecturas limpias y ordena el BAM resultante.
rule bwa_mem_sort:
    input:
        r1="results/trimmed/{sample}_R1.fastq.gz",
        r2="results/trimmed/{sample}_R2.fastq.gz",
        ref=config["ref_genome"],
        # Dependencia explícita en la regla de indexación.
        indexed=config["ref_genome"] + ".bwa_indexed"
    output:
        bam="results/aligned/{sample}.sorted.bam"
    params:
        # Añadir 'read groups' es una buena práctica y requerido por GATK.
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

# ==============================================================================
# PASO 4: PROCESAMIENTO POST-ALINEAMIENTO
# ==============================================================================
# Marca duplicados de PCR con GATK.
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

# Indexa el BAM final (deduplicado) para acceso rápido.
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

# ==============================================================================
# PASO 5: LLAMADA DE VARIANTES
# ==============================================================================
# Llama variantes usando GATK HaplotypeCaller.
rule haplotype_caller:
    input:
        bam="results/aligned/{sample}.dedup.bam",
        bai="results/aligned/{sample}.dedup.bam.bai",
        ref=config["ref_genome"]
    output:
        vcf="results/variants/{sample}.vcf.gz"
    params:
        # Configura la memoria RAM para la JVM de GATK.
        # Usa el valor pasado desde la línea de comando (por run_pipeline.sh).
        # Si no se pasa, usa un valor seguro por defecto de 4GB.
        java_opts=f"-Xmx{config.get('gatk_ram_gb', 4)}g"
    log:
        "results/logs/haplotype_caller/{sample}.log"
    container:
        "docker://broadinstitute/gatk:4.2.6.1"
    shell:
        "gatk --java-options '{params.java_opts}' HaplotypeCaller "
        "-R {input.ref} "
        "-I {input.bam} "
        "-O {output.vcf} > {log} 2>&1"

# ==============================================================================
# PASO 6: ANOTACIÓN DE VARIANTES
# ==============================================================================
# Anota el VCF con SnpEff para predecir el impacto de las variantes.
rule snpeff_annotate:
    input:
        vcf="results/variants/{sample}.vcf.gz",
        ref=config["ref_genome"]
    output:
        vcf="results/annotated/{sample}.ann.vcf",
        html="results/annotated/{sample}.snpeff.html"
    params:
        # El nombre de la base de datos se lee desde el config.yaml.
        db=config["snpeff"]["database"]
    log:
        "results/logs/snpeff/{sample}.log"
    container:
        "docker://snpeff/snpeff:5.1d-1"
    shell:
        "snpEff ann -v {params.db} {input.vcf} > {output.vcf} "
        "-stats {output.html} 2> {log}"
