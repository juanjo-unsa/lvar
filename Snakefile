# ==============================================================================
# Snakefile para LVAR (Versión Final y Robusta con SnpEff)
# ==============================================================================

import pandas as pd

# --- Carga de Configuracion ---
configfile: "config/config.yaml"
SAMPLES = pd.read_csv(config['samples'], sep='\t').set_index('sample', drop=False)

# --- Rutas y Parametros Globales (dentro del contenedor) ---
REF_GENOME = "/opt/db/genome.fasta"
SNPEFF_DB = "LbraziliensisMHOMBR75M2904"

# ==============================================================================
# REGLA FINAL (TARGET)
# ==============================================================================
rule all:
    input:
        expand("results/annotated/{sample}.ann.vcf", sample=SAMPLES.index),
        "results/qc/multiqc_report.html"

# ... (Las reglas fastqc, multiqc, fastp, bwa_index, bwa_mem_sort, mark_duplicates, samtools_index, haplotype_caller no cambian, solo asegúrate de que usan la variable REF_GENOME) ...

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
