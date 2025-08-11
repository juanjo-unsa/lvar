# ==============================================================================
# Dockerfile para LVAR (Versión Definitiva con Snakemake completo)
# ==============================================================================

# Usar una base de Debian para un mejor control de las dependencias de Perl
FROM debian:bullseye-slim

# --- Metadata ---
LABEL maintainer="Juanjo <tu.email@ejemplo.com>"
LABEL description="Entorno completo para el pipeline LVAR con Ensembl VEP."

# --- Instalar dependencias del sistema ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    unzip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libdbi-perl \
    libdbd-mysql-perl \
    perl \
    cpanminus \
    procps \
    ca-certificates && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# --- Instalar Mambaforge ---
RUN wget --quiet --no-check-certificate https://github.com/conda-forge/miniforge/releases/download/24.1.2-0/Mambaforge-24.1.2-0-Linux-x86_64.sh -O mambaforge.sh && \
    bash mambaforge.sh -b -p /opt/conda && \
    rm mambaforge.sh
ENV PATH /opt/conda/bin:$PATH

# --- Instalar herramientas bioinformáticas con Mamba (CORRECCIÓN APLICADA AQUÍ) ---
# Se instala el paquete 'snakemake' completo en lugar de 'snakemake-minimal'
# para asegurar que dependencias como 'pandas' estén incluidas.
RUN mamba install -n base -c conda-forge -c bioconda -y \
    snakemake \
    fastqc \
    multiqc \
    fastp \
    bwa \
    samtools \
    gatk4 \
    htslib \
    ensembl-vep && \
    mamba clean --all -y

# --- Configurar Ensembl VEP con un caché autocontenido ---
# (Esta sección no necesita cambios)
ENV VEP_CACHE_DIR /opt/vep_cache
ENV VEP_SPECIES leishmania_braziliensis_mhomb_br_75_m2904
ENV VEP_ASSEMBLY ASM244v1
ENV VEP_VERSION 112
ENV VEP_FULL_VERSION ${VEP_VERSION}_${VEP_ASSEMBLY}

RUN mkdir -p ${VEP_CACHE_DIR} && \
    wget -P ${VEP_CACHE_DIR} https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-58/fasta/leishmania_braziliensis_mhomb_br_75_m2904/dna/Leishmania_braziliensis_mhomb_br_75_m2904.ASM244v1.dna.toplevel.fa.gz && \
    wget -P ${VEP_CACHE_DIR} https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-58/gff3/leishmania_braziliensis_mhomb_br_75_m2904/Leishmania_braziliensis_mhomb_br_75_m2904.ASM244v1.58.gff3.gz && \
    gunzip ${VEP_CACHE_DIR}/*.gz && \
    vep_install --AUTO c --SPECIES ${VEP_SPECIES} --ASSEMBLY ${VEP_ASSEMBLY} \
        --CACHEDIR ${VEP_CACHE_DIR} --NO_UPDATE --QUIET && \
    samtools faidx ${VEP_CACHE_DIR}/${VEP_SPECIES}/${VEP_FULL_VERSION}/Leishmania_braziliensis_mhomb_br_75_m2904.ASM244v1.dna.toplevel.fa

# --- Configurar el entorno de trabajo ---
WORKDIR /pipeline

# --- Definir el punto de entrada ---
ENTRYPOINT ["snakemake"]
CMD ["--help"]
