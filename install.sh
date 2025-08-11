#!/bin/bash
# -----------------------------------------------------------------------------
# install.sh - Script de configuración interactivo para el proyecto 'lvar'
#
# Este script configura el entorno, genera los archivos de configuración
# y prepara el proyecto para su ejecución.
# -----------------------------------------------------------------------------

# --- Colores y Funciones ---
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
echo -e "${GREEN}--- Configurando el entorno del proyecto 'lvar' ---${NC}"

# --- PASO 0: Verificar ubicación y dependencias ---
if [ ! -f "Snakefile" ]; then echo -e "${RED}[ERROR] 'Snakefile' no encontrado. Ejecuta este script desde la raíz del proyecto 'lvar'.${NC}"; exit 1; fi
echo -e "\n${GREEN}1. Verificando dependencias del sistema...${NC}"
check_dependency() { local cmd=$1; local msg=$2; echo -n "   Verificando '$cmd'... "; if ! command -v "$cmd" &> /dev/null; then echo -e "${RED}NO ENCONTRADO\n${YELLOW}Por favor, instala '$cmd'. ${msg}${NC}"; exit 1; fi; echo -e "${GREEN}OK${NC}"; }
check_dependency "docker" "Visita https://docs.docker.com/get-docker/"
check_dependency "snakemake" "Recomendado vía Conda: 'conda install -c bioconda snakemake'"

# --- PASO 2: Crear estructura de directorios ---
echo -e "\n${GREEN}2. Creando la estructura de directorios del proyecto...${NC}"
mkdir -p config data/{raw_fastq,reference} results/{qc,trimmed,aligned,variants,annotated,logs} docs
echo "   [OK] Estructura de directorios creada."

# --- PASO 3: Generar config/samples.tsv ---
echo -e "\n${GREEN}3. Configurando las muestras del análisis...${NC}"
declare -a SAMPLES_ARRAY
echo -e "${YELLOW}Ahora ingresa los nombres de tus muestras, uno por uno.${NC}"
echo "El 'nombre de la muestra' es el prefijo de tus archivos FASTQ (ej: para 'susceptible_R1.fastq.gz', el nombre es 'susceptible')."
while true; do
    read -p "Ingresa un nombre de muestra (o presiona Enter para terminar): " sample_name
    if [ -z "$sample_name" ]; then
        break
    fi
    SAMPLES_ARRAY+=("$sample_name")
done

if [ ${#SAMPLES_ARRAY[@]} -eq 0 ]; then
    echo -e "${RED}[ERROR] No se ingresó ninguna muestra. El pipeline no puede continuar.${NC}"
    exit 1
fi

echo "sample" > config/samples.tsv
for sample in "${SAMPLES_ARRAY[@]}"; do
    echo "$sample" >> config/samples.tsv
done
echo "   [OK] Archivo 'config/samples.tsv' generado con ${#SAMPLES_ARRAY[@]} muestras."

# --- PASO 4: Generar config/config.yaml ---
echo -e "\n${GREEN}4. Configurando los parámetros del pipeline...${NC}"
read -p "Ingresa el nombre del archivo del genoma de referencia [L_braziliensis_ref.fasta]: " REF_FILENAME
REF_FILENAME=${REF_FILENAME:-"L_braziliensis_ref.fasta"}

echo -e "${YELLOW}Ingresa el nombre de la base de datos de SnpEff para tu organismo.${NC}"
echo "Este nombre debe coincidir con una base de datos instalada en SnpEff."
echo "Ejemplo para L. braziliensis MHOM/BR/75/M2904: LbraziliensisMHOMBR75M2904"
read -p "Nombre de la base de datos de SnpEff: " SNPEFF_DB

if [ -z "$SNPEFF_DB" ]; then
    echo -e "${RED}[ERROR] El nombre de la base de datos de SnpEff es obligatorio.${NC}"; exit 1
fi

cat > config/config.yaml <<- EOM
# Archivo de configuración generado por install.sh

# Ruta al archivo de muestras (generado automáticamente)
samples: "config/samples.tsv"

# Ruta al genoma de referencia (relativa a la raíz del proyecto)
ref_genome: "data/reference/${REF_FILENAME}"

# Parámetros de SnpEff
snpeff:
  database: "${SNPEFF_DB}"
EOM
echo "   [OK] Archivo 'config/config.yaml' generado."

# --- PASO 5: Configurar recursos y generar script de ejecución ---
echo -e "\n${GREEN}5. Configurando los recursos del sistema...${NC}"
# (Esta sección no cambia)
CORES_TOTAL=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 8)
RAM_SUGGESTED=$(( $(free -g | awk '/^Mem:/{print $2}' 2>/dev/null || echo 16) * 80 / 100 ))
read -p "Número de cores/threads a usar [Sugerido: $CORES_TOTAL]: " USER_CORES
read -p "Memoria RAM máxima para GATK (en GB) [Sugerido: ${RAM_SUGGESTED}]: " USER_RAM
USER_CORES=${USER_CORES:-$CORES_TOTAL}; USER_RAM=${USER_RAM:-$RAM_SUGGESTED}
if ! [[ "$USER_CORES" =~ ^[0-9]+$ && "$USER_RAM" =~ ^[0-9]+$ ]]; then echo -e "${RED}Cores y RAM deben ser números.${NC}"; exit 1; fi
cat > run_pipeline.sh <<- EOM
#!/bin/bash
echo "Iniciando pipeline 'lvar' con ${USER_CORES} cores y GATK con ${USER_RAM}GB RAM..."
snakemake --use-software-stack --cores \${USER_CORES} --config gatk_ram_gb=\${USER_RAM} "\$@"
EOM
chmod +x run_pipeline.sh
echo "   [OK] Script de ejecución './run_pipeline.sh' creado."

# --- PASO 6: Copiar y verificar archivos de datos ---
echo -e "\n${GREEN}6. Copiando y verificando los archivos de datos...${NC}"
read -e -p "Proporciona la ruta a la carpeta donde tienes tus archivos de entrada: " SOURCE_DIR
if [ -d "$SOURCE_DIR" ]; then
    echo "   Copiando genoma de referencia..."
    if [ -f "$SOURCE_DIR/$REF_FILENAME" ]; then cp "$SOURCE_DIR/$REF_FILENAME" "data/reference/"; fi
    echo "   Copiando lecturas de las muestras..."
    for sample in "${SAMPLES_ARRAY[@]}"; do
        for read_pair in R1 R2; do
            fastq_file="${sample}_${read_pair}.fastq.gz"
            if [ -f "$SOURCE_DIR/$fastq_file" ]; then cp "$SOURCE_DIR/$fastq_file" "data/raw_fastq/"; fi
        done
    done
else
    echo -e "${RED}[ADVERTENCIA] El directorio '$SOURCE_DIR' no existe. Deberás mover los archivos manualmente.${NC}"
fi

echo -e "\n${GREEN}Verificación final:${NC}"
echo "-------------------------------------------------------------------------"
check_file() {
    local file_path=$1; local missing_msg=$2
    while [ ! -f "$file_path" ]; do echo -e "${RED}[ACCIÓN] Falta archivo: $file_path${NC}\n${YELLOW}$missing_msg${NC}"; read -p "Presiona [Enter] tras colocar el archivo..."; done
    echo -e "${GREEN}[OK] Archivo verificado: $file_path${NC}"
}
check_file "data/reference/$REF_FILENAME" "Mueve el genoma de referencia a 'data/reference/'."
for sample in "${SAMPLES_ARRAY[@]}"; do
    check_file "data/raw_fastq/${sample}_R1.fastq.gz" "Mueve las lecturas R1 de la muestra '$sample' a 'data/raw_fastq/'."
    check_file "data/raw_fastq/${sample}_R2.fastq.gz" "Mueve las lecturas R2 de la muestra '$sample' a 'data/raw_fastq/'."
done
echo "-------------------------------------------------------------------------"

# --- Conclusión ---
echo -e "\n${GREEN}¡CONFIGURACIÓN COMPLETADA!${NC}"
echo "El entorno del proyecto 'lvar' está listo para ser utilizado."
echo "Para iniciar el análisis, ejecuta:"
echo -e "${YELLOW}./run_pipeline.sh${NC}"
