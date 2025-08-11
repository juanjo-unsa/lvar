# ==============================================================================
# Dockerfile para LVAR (Versión Definitiva y Robusta con SnpEff)
# ==============================================================================

# Usar Mambaforge para evitar problemas de ToS de Conda
FROM condaforge/mambaforge:latest

# --- Metadata ---
LABEL maintainer="Juanjo <tu.email@ejemplo.com>"
LABEL description="Entorno completo para el pipeline LVAR con SnpEff."

# --- Instalar herramientas bioinformáticas con Mamba ---
RUN mamba install -n base -c conda-forge -c bioconda -y \
    snakemake \
    fastqc \
    multiqc \
    fastp \
    bwa \
    samtools \
    gatk4 \
    snpeff \
    wget && \
    mamba clean --all -y

# --- Configurar SnpEff y el Genoma de Referencia ---
# Usamos la base de datos de referencia principal para L. braziliensis.
ENV SNPEFF_DB LbraziliensisMHOMBR75M2904
# Ruta donde guardaremos el genoma dentro del contenedor
ENV GENOME_PATH "/opt/db/genome.fasta"

RUN mkdir -p /opt/db && \
    # Descargar la base de datos de SnpEff. El flag -v fuerza la búsqueda web.
    snpEff download -v ${SNPEFF_DB} && \
    # Extraer la ruta al FASTA desde la configuración de SnpEff para asegurar coherencia
    # y descargarlo a nuestra ubicación estándar.
    GENOME_URL=$(grep "${SNPEFF_DB}.genome" $(find /opt/conda/share -name snpEff.config) | cut -d'=' -f2 | tr -d '[:space:]') && \
    wget -O ${GENOME_PATH}.gz "https://initial-pre-sept-2023.s3.eu-west-1.amazonaws.com/downloads.sourceforge.net/project/snpeff/databases/v5_1/${GENOME_URL}" && \
    gunzip ${GENOME_PATH}.gz

# --- Configurar el entorno de trabajo ---
WORKDIR /pipeline

# --- Definir el punto de entrada ---
ENTRYPOINT ["snakemake"]
CMD ["--help"]
