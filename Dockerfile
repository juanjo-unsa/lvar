# ==============================================================================
# Dockerfile para el proyecto LVAR (Versión Final)
# ==============================================================================

# Usar una base con Miniconda para gestionar los paquetes
FROM continuumio/miniconda3:latest

# --- Metadata ---
LABEL maintainer="Juan José Aguirre"
LABEL description="Entorno completo para el pipeline LVAR de análisis de variantes."

# --- Argumento de Construcción para SnpEff ---
# Esto permite pasar el nombre de la base de datos durante la construcción.
# Se establece un valor por defecto por si no se pasa ninguno.
ARG SNPEFF_DB=LbraziliensisMHOMBR75M2904

# --- Instalar dependencias del sistema y Mamba ---
RUN apt-get update && apt-get install -y procps && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN conda install -n base -c conda-forge -y mamba

# --- Instalar herramientas bioinformáticas con Mamba ---
RUN mamba install -n base -c conda-forge -c bioconda -y \
    snakemake-minimal \
    fastqc \
    multiqc \
    fastp \
    bwa \
    samtools \
    gatk4 \
    snpeff && \
    mamba clean --all -y

# --- Configurar SnpEff ---
# Descarga la base de datos especificada por el argumento de construcción.
RUN echo "Descargando la base de datos de SnpEff: ${SNPEFF_DB}" && \
    snpEff download ${SNPEFF_DB}

# --- Configurar el entorno de trabajo ---
WORKDIR /pipeline

# --- Definir el punto de entrada ---
ENTRYPOINT ["snakemake"]
CMD ["--help"]
