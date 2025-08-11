# ==============================================================================
# Dockerfile para el proyecto LVAR (Versión Final y Robusta)
#
# Construye una imagen con Ensembl VEP y un caché pre-configurado para
# Leishmania braziliensis, asegurando la máxima reproducibilidad.
# ==============================================================================

# Usar una base de Debian para un mejor control de las dependencias de Perl
FROM debian:bullseye-slim

# --- Metadata ---
LABEL maintainer="Juanjo <tu.email@ejemplo.com>"
LABEL description="Entorno completo para el pipeline LVAR con Ensembl VEP."

# --- Instalar dependencias del sistema ---
# VEP necesita Perl, git, y otras herramientas de compilación y librerías.
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

# --- Instalar Miniconda y Mamba ---
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh
ENV PATH /opt/conda/bin:$PATH
RUN conda install -n base -c conda-forge -y mamba

# --- Instalar herramientas bioinformáticas ---
RUN mamba install -n base -c conda-forge -c bioconda -y \
    snakemake-minimal \
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
# Este es el paso más largo, pero crucial. Descarga y prepara todo lo necesario.
# El nombre de la especie es genérico, pero los archivos son específicos.
# Leishmania braziliensis (MHOM/BR/75/M2904) de Ensembl Protists release 58
ENV VEP_CACHE_DIR /opt/vep_cache
ENV VEP_SPECIES leishmania_braziliensis_mhomb_br_75_m2904
ENV VEP_ASSEMBLY ASM244v1
ENV VEP_VERSION 112
ENV VEP_FULL_VERSION ${VEP_VERSION}_${VEP_ASSEMBLY}

RUN mkdir -p ${VEP_CACHE_DIR} && \
    # Descargar el genoma (FASTA) y la anotación (GFF3)
    wget -P ${VEP_CACHE_DIR} https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-58/fasta/leishmania_braziliensis_mhomb_br_75_m2904/dna/Leishmania_braziliensis_mhomb_br_75_m2904.ASM244v1.dna.toplevel.fa.gz && \
    wget -P ${VEP_CACHE_DIR} https://ftp.ensemblgenomes.ebi.ac.uk/pub/protists/release-58/gff3/leishmania_braziliensis_mhomb_br_75_m2904/Leishmania_braziliensis_mhomb_br_75_m2904.ASM244v1.58.gff3.gz && \
    # Descomprimir y preparar
    gunzip ${VEP_CACHE_DIR}/*.gz && \
    # Construir el caché de VEP
    vep_install --AUTO c --SPECIES ${VEP_SPECIES} --ASSEMBLY ${VEP_ASSEMBLY} \
        --CACHEDIR ${VEP_CACHE_DIR} --NO_UPDATE --QUIET && \
    # Indexar el FASTA para samtools
    samtools faidx ${VEP_CACHE_DIR}/leishmania_braziliensis_mhomb_br_75_m2904/112_ASM244v1/Leishmania_braziliensis_mhomb_br_75_m2904.ASM244v1.dna.toplevel.fa

# --- Configurar el entorno de trabajo ---
WORKDIR /pipeline

# --- Definir el punto de entrada ---
ENTRYPOINT ["snakemake"]
CMD ["--help"]
