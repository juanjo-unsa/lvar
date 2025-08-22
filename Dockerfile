# ==============================================================================
# Dockerfile para el pipeline LVAR (FINAL Y A PRUEBA DE FALLOS)
# ==============================================================================

# Usar una base limpia de Debian (stable) para asegurar disponibilidad de utilidades basicas.
FROM debian:stable

# --- ARGUMENTOS DE CONSTRUCCION ---
# Estos argumentos seran pasados por el script setup.sh
ARG FASTA_URL_ARG
ARG GFF_URL_ARG

# --- Metadata ---
LABEL maintainer="Juanjo <tu.email@ejemplo.com>"
LABEL description="Entorno completo para el pipeline LVAR con SnpEff y bcftools."

# --- Instalar dependencias del sistema, incluyendo `less`. ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    build-essential \
    procps \
    less \
    ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# --- Instalar Mambaforge (gestor de paquetes) ---
RUN wget --quiet --no-check-certificate https://github.com/conda-forge/miniforge/releases/download/24.1.2-0/Mambaforge-24.1.2-0-Linux-x86_64.sh -O mambaforge.sh && \
    bash mambaforge.sh -b -p /opt/conda && \
    rm mambaforge.sh
ENV PATH /opt/conda/bin:$PATH

# --- Instalar todas las herramientas bioinformaticas. ---
RUN mamba install -n base -c conda-forge -c bioconda -y \
    snakemake \
    fastqc \
    multiqc \
    fastp \
    bwa \
    samtools \
    bcftools \
    gatk4 \
    snpeff && \
    mamba clean --all -y

# --- Construir manualmente la base de datos de SnpEff. ---
RUN \
    # Definir variables internas usando los argumentos de construccion
    DB_NAME="Lbraziliensis_custom_db" && \
    FASTA_URL="${FASTA_URL_ARG}" && \
    GFF_URL="${GFF_URL_ARG}" && \
    GENOME_PATH="/opt/db/genome.fasta" && \
    # Encontrar la ruta de SnpEff dinamicamente
    SNPEFF_CONFIG_FILE=$(find /opt/conda/share -name snpEff.config) && \
    SNPEFF_DIR=$(dirname "${SNPEFF_CONFIG_FILE}") && \
    SNPEFF_DATA_DIR="${SNPEFF_DIR}/data" && \
    \
    # Crear directorios
    mkdir -p "${SNPEFF_DATA_DIR}/${DB_NAME}" && \
    mkdir -p /opt/db && \
    \
    # Descargar los archivos de referencia de TriTrypDB
    wget -O "${SNPEFF_DATA_DIR}/${DB_NAME}/genes.gff" "${GFF_URL}" && \
    wget -O "${GENOME_PATH}" "${FASTA_URL}" && \
    \
    # Preparar archivos para SnpEff
    cp "${GENOME_PATH}" "${SNPEFF_DATA_DIR}/${DB_NAME}/sequences.fa" && \
    \
    # Anadir nuestra base de datos personalizada a la configuracion de SnpEff
    echo "" >> "${SNPEFF_CONFIG_FILE}" && \
    echo "# Base de datos manual de L. braziliensis (proporcionada por el usuario)" >> "${SNPEFF_CONFIG_FILE}" && \
    echo "${DB_NAME}.genome : Leishmania braziliensis (custom)" >> "${SNPEFF_CONFIG_FILE}" && \
    \
    # Construir la base de datos de SnpEff
    snpEff build -gff3 -v "${DB_NAME}" -noCheckCds -noCheckProtein

# Configurar el entorno de trabajo y el punto de entrada para Snakemake.
WORKDIR /pipeline
ENTRYPOINT ["snakemake"]
CMD ["--help"]
