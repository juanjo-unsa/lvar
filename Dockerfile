# ==============================================================================
# Dockerfile para LVAR (Versión Final Definitiva - Flags SnpEff Corregidos)
# ==============================================================================

# Usar Mambaforge para evitar problemas de ToS de Conda y tener Mamba por defecto
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

# --- Configuración Manual de SnpEff y el Genoma de Referencia ---
# Descargamos el Genoma y el GFF, y construimos la DB con los flags correctos
# para omitir las comprobaciones de sanidad.

# Variables para las URLs (VERIFICADAS) y nombres
ENV DB_NAME "Lbraziliensis_2019_manual"
ENV FASTA_URL "https://tritrypdb.org/common/downloads/Current_Release/LbraziliensisMHOMBR75M2904_2019/fasta/data/TriTrypDB-68_LbraziliensisMHOMBR75M2904_2019_Genome.fasta"
ENV GFF_URL "https://tritrypdb.org/common/downloads/Current_Release/LbraziliensisMHOMBR75M2904_2019/gff/data/TriTrypDB-68_LbraziliensisMHOMBR75M2904_2019.gff"

# Ruta al directorio de datos de SnpEff dentro del contenedor
ENV SNPEFF_DATA_DIR "/opt/conda/share/snpeff-5.2-1/data"
# Ruta donde guardaremos nuestro genoma para el alineamiento
ENV GENOME_PATH "/opt/db/genome.fasta"

RUN mkdir -p ${SNPEFF_DATA_DIR}/${DB_NAME} && \
    mkdir -p /opt/db && \
    # Descargar el genoma y la anotación
    wget -O ${SNPEFF_DATA_DIR}/${DB_NAME}/genes.gff "${GFF_URL}" && \
    wget -O ${GENOME_PATH} "${FASTA_URL}" && \
    # Copiar el genoma a la carpeta de SnpEff con el nombre que espera
    cp ${GENOME_PATH} ${SNPEFF_DATA_DIR}/${DB_NAME}/sequences.fa && \
    # Añadir la configuración de nuestra base de datos manual al snpEff.config
    echo "" >> /opt/conda/share/snpeff-5.2-1/snpEff.config && \
    echo "# Manual database for L. braziliensis 2019 from TriTrypDB" >> /opt/conda/share/snpeff-5.2-1/snpEff.config && \
    echo "${DB_NAME}.genome : Leishmania braziliensis 2019 (manual)" >> /opt/conda/share/snpeff-5.2-1/snpEff.config && \
    # Construir la base de datos, usando los flags correctos para esta versión de SnpEff
    snpEff build -gff3 -v ${DB_NAME} -noCheckCds -noCheckProtein

# --- Configurar el entorno de trabajo ---
WORKDIR /pipeline

# --- Definir el punto de entrada ---
ENTRYPOINT ["snakemake"]
CMD ["--help"]
