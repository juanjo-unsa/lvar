# ==============================================================================
# Dockerfile para LVAR (Versión Final - Construcción Completa y Verificada de SnpEff)
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
# Este es el método más robusto: descargamos TODOS los archivos y construimos la DB.

# Variables para las URLs (PROPORCIONADAS Y VERIFICADAS) y nombres
ENV DB_NAME "Lbraziliensis_2019_manual"
ENV BASE_URL "https://tritrypdb.org/common/downloads/Current_Release/LbraziliensisMHOMBR75M2904_2019"
ENV GFF_URL "${BASE_URL}/gff/data/TriTrypDB-68_LbraziliensisMHOMBR75M2904_2019.gff"
ENV FASTA_URL "${BASE_URL}/fasta/data/TriTrypDB-68_LbraziliensisMHOMBR75M2904_2019_Genome.fasta"
ENV CDS_URL "${BASE_URL}/fasta/data/TriTrypDB-68_LbraziliensisMHOMBR75M2904_2019_AnnotatedCDSs.fasta"
ENV PROTEIN_URL "${BASE_URL}/fasta/data/TriTrypDB-68_LbraziliensisMHOMBR75M2904_2019_AnnotatedProteins.fasta"

# Ruta al directorio de datos de SnpEff dentro del contenedor
ENV SNPEFF_DATA_DIR "/opt/conda/share/snpeff-5.2-1/data"
# Ruta donde guardaremos nuestro genoma para el alineamiento
ENV GENOME_PATH "/opt/db/genome.fasta"

RUN mkdir -p ${SNPEFF_DATA_DIR}/${DB_NAME} && \
    mkdir -p /opt/db && \
    # Descargar los 4 archivos necesarios
    wget -O ${SNPEFF_DATA_DIR}/${DB_NAME}/genes.gff "${GFF_URL}" && \
    wget -O ${GENOME_PATH} "${FASTA_URL}" && \
    wget -O ${SNPEFF_DATA_DIR}/${DB_NAME}/cds.fa "${CDS_URL}" && \
    wget -O ${SNPEFF_DATA_DIR}/${DB_NAME}/protein.fa "${PROTEIN_URL}" && \
    # Copiar el genoma a la carpeta de SnpEff con el nombre que espera
    cp ${GENOME_PATH} ${SNPEFF_DATA_DIR}/${DB_NAME}/sequences.fa && \
    # Añadir la configuración de nuestra base de datos manual al snpEff.config
    echo "" >> /opt/conda/share/snpeff-5.2-1/snpEff.config && \
    echo "# Manual database for L. braziliensis 2019 from TriTrypDB" >> /opt/conda/share/snpeff-5.2-1/snpEff.config && \
    echo "${DB_NAME}.genome : Leishmania braziliensis 2019 (manual)" >> /opt/conda/share/snpeff-5.2-1/snpEff.config && \
    # Construir la base de datos binaria de SnpEff. Ya no se necesita -noCheck.
    snpEff build -gff3 -v ${DB_NAME}

# --- Configurar el entorno de trabajo ---
WORKDIR /pipeline

# --- Definir el punto de entrada ---
ENTRYPOINT ["snakemake"]
CMD ["--help"]
