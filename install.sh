#!/bin/bash
# =============================================================================
# LVAR - Universal Installer Script (Versión Final)
# =============================================================================

# --- Configuración y Colores ---
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
PROJECT_NAME="lvar_project"
REPO_URL="https://github.com/juanjo-unsa/lvar.git"
DOCKER_IMAGE_TAG="lvar-vep-pipeline:latest"

echo -e "${GREEN}--- Iniciando la Instalación del Pipeline LVAR (con Ensembl VEP) ---${NC}"

# --- PASO 1: Verificar Dependencias ---
echo -e "\n${GREEN}1. Verificando dependencias del sistema...${NC}"
check_dependency() {
    local cmd=$1; local msg=$2; echo -n "   Verificando '$cmd'... ";
    if ! command -v "$cmd" &> /dev/null; then echo -e "${RED}NO ENCONTRADO\n${YELLOW}Por favor, instala '$cmd'. ${msg}${NC}"; exit 1; fi;
    echo -e "${GREEN}OK${NC}";
}
check_dependency "git" "Ej: 'sudo apt-get install git'"
check_dependency "docker" "Visita https://docs.docker.com/get-docker/"

# --- PASO 2: Clonar Repositorio y Crear Estructura ---
echo -e "\n${GREEN}2. Descargando el código fuente del pipeline...${NC}"
if [ -d "$PROJECT_NAME" ]; then echo -e "${RED}[ERROR] El directorio '$PROJECT_NAME' ya existe.${NC}"; exit 1; fi
git clone "$REPO_URL" "$PROJECT_NAME" || { echo -e "${RED}Falló la clonación.${NC}"; exit 1; }
cd "$PROJECT_NAME" || exit
mkdir -p config data/raw_fastq results/{qc,trimmed,aligned,variants,annotated,logs}
echo "   [OK] Proyecto clonado y estructura de carpetas creada en: $(pwd)"

# --- PASO 3: Construir la Imagen de Docker ---
echo -e "\n${GREEN}3. Construyendo la imagen de Docker del pipeline...${NC}"
echo -e "${YELLOW}Este es el paso más largo y puede tardar más de 30 minutos.${NC}"
echo -e "${YELLOW}Descargará e instalará todas las herramientas y el genoma de referencia.${NC}"
if docker build -t "$DOCKER_IMAGE_TAG" .; then
    echo -e "   ${GREEN}[OK] Imagen Docker '$DOCKER_IMAGE_TAG' construida con éxito.${NC}"
else
    echo -e "${RED}[ERROR] Falló la construcción de la imagen de Docker. Revisa los mensajes de error.${NC}"; exit 1
fi

# --- PASO 4: Configuración Interactiva ---
echo -e "\n${GREEN}4. Configurando las muestras y parámetros del pipeline...${NC}"
declare -a SAMPLES_ARRAY
echo -e "${YELLOW}Ingresa los nombres de tus muestras (prefijos de los archivos FASTQ). Presiona Enter para terminar.${NC}"
while true; do read -p "Nombre de muestra: " sample_name; [ -z "$sample_name" ] && break; SAMPLES_ARRAY+=("$sample_name"); done
if [ ${#SAMPLES_ARRAY[@]} -eq 0 ]; then echo -e "${RED}No se ingresaron muestras. Abortando.${NC}"; exit 1; fi

echo "sample" > config/samples.tsv; for s in "${SAMPLES_ARRAY[@]}"; do echo "$s" >> config/samples.tsv; done
# El config.yaml ahora es mínimo, ya que todo está en el contenedor.
cat > config/config.yaml <<- EOM
# Ruta al archivo de muestras (generado automáticamente)
samples: "config/samples.tsv"
EOM
echo "   [OK] Archivos 'config/samples.tsv' y 'config/config.yaml' generados."

# --- PASO 5: Generar Script de Ejecución ---
echo -e "\n${GREEN}5. Creando script de ejecución personalizado...${NC}"
CORES_TOTAL=$(nproc 2>/dev/null || echo 8)
RAM_SUGGESTED=$(( $(free -g | awk '/^Mem:/{print $2}' 2>/dev/null || echo 16) * 80 / 100 ))
read -p "Número de cores a usar [Sugerido: $CORES_TOTAL]: " USER_CORES; USER_CORES=${USER_CORES:-$CORES_TOTAL}
read -p "RAM para GATK (GB) [Sugerido: ${RAM_SUGGESTED}]: " USER_RAM; USER_RAM=${USER_RAM:-$RAM_SUGGESTED}
cat > run_docker.sh <<- EOM
#!/bin/bash
# Lanza el pipeline dentro del contenedor Docker pre-construido.

echo "Iniciando pipeline LVAR con ${USER_CORES} cores y GATK con ${USER_RAM}GB RAM..."
docker run --rm -it \\
    -v "\$(pwd)/config:/pipeline/config" \\
    -v "\$(pwd)/data:/pipeline/data" \\
    -v "\$(pwd)/results:/pipeline/results" \\
    -v "\$(pwd)/Snakefile:/pipeline/Snakefile" \\
    "${DOCKER_IMAGE_TAG}" \\
    --cores \${USER_CORES} \\
    --config gatk_ram_gb=\${USER_RAM} \\
    "\$@"
EOM
chmod +x run_docker.sh
echo "   [OK] Script de ejecución './run_docker.sh' creado."

# --- PASO 6: Organizar Archivos de Datos ---
echo -e "\n${GREEN}6. Organizando los archivos de datos FASTQ...${NC}"
read -e -p "Proporciona la RUTA ABSOLUTA a tus datos de origen (la carpeta con los .fastq.gz): " SOURCE_DIR
if [ -d "$SOURCE_DIR" ]; then
    for s in "${SAMPLES_ARRAY[@]}"; do
        cp "$SOURCE_DIR/${s}_R1.fastq.gz" "data/raw_fastq/" 2>/dev/null || echo -e "${RED}No se pudo copiar ${s}_R1.fastq.gz${NC}"
        cp "$SOURCE_DIR/${s}_R2.fastq.gz" "data/raw_fastq/" 2>/dev/null || echo -e "${RED}No se pudo copiar ${s}_R2.fastq.gz${NC}"
    done
    echo "   [OK] Intento de copia finalizado. Por favor, verifica los archivos en 'data/raw_fastq/'."
else
    echo -e "${RED}Directorio no encontrado. Deberás mover los archivos FASTQ manualmente a data/raw_fastq/.${NC}"
fi

# --- Conclusión ---
echo -e "\n${GREEN}¡INSTALACIÓN COMPLETA!${NC}"
echo "El proyecto LVAR ha sido instalado en: ${YELLOW}$(pwd)${NC}"
echo "Para ejecutar el análisis, usa el script:"
echo -e "${GREEN}./run_docker.sh${NC}"
