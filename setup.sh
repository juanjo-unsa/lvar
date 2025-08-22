#!/bin/bash
# =============================================================================
# LVAR Project Setup Script (setup.sh) - Versión Final Interactiva
# =============================================================================

# --- Configuración y Colores ---
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
DOCKER_IMAGE_TAG="lvar-snpeff-pipeline:latest"

# --- Función para obtener y validar una URL del usuario (sin valor por defecto) ---
get_valid_url() {
    local prompt_message=$1
    local user_url

    while true; do
        read -e -p "$(echo -e "${YELLOW}${prompt_message}${NC} ")" user_url
        if [ -z "$user_url" ]; then
            echo -e "${RED}[ERROR] El enlace no puede estar vacío. Por favor, inténtalo de nuevo.${NC}"
            continue
        fi

        echo -e "${GREEN}--> Verificando enlace...${NC}"
        if wget --spider --quiet "$user_url"; then
            echo -e "   ${GREEN}[OK] El enlace es válido.${NC}"
            echo "$user_url"
            break
        else
            echo -e "   ${RED}[ERROR] El enlace no es válido o el servidor no responde. Por favor, inténtalo de nuevo.${NC}"
        fi
    done
}

echo -e "${GREEN}--- Iniciando la Configuración del Proyecto LVAR ---${NC}"

# --- PASO 1: Verificar Requisitos ---
if [ ! -f "Dockerfile" ] || [ ! -f "Snakefile" ]; then echo -e "${RED}ERROR: 'Dockerfile' o 'Snakefile' no encontrados. Ejecuta desde la raíz del repo.${NC}"; exit 1; fi
echo -e "\n${GREEN}1. Verificando la presencia de Docker...${NC}"
if ! command -v docker &> /dev/null; then echo -e "${RED}ERROR: Docker no instalado.${NC}"; exit 1; fi
echo "   [OK] Docker está instalado."

# --- PASO 2: Configurar URLs de Referencia ---
echo -e "\n${GREEN}2. Configurando enlaces para el genoma y la anotación...${NC}"
FASTA_URL=$(get_valid_url "Por favor, introduce la URL para el genoma (.fasta):")
GFF_URL=$(get_valid_url "Ahora, introduce la URL para la anotación (.gff):")

# --- PASO 3: Construir Imagen Docker ---
echo -e "\n${GREEN}3. Verificando y construyendo la imagen de Docker...${NC}"
if [[ "$(docker images -q ${DOCKER_IMAGE_TAG} 2> /dev/null)" == "" ]]; then
    echo -e "${YELLOW}La imagen '${DOCKER_IMAGE_TAG}' no existe. Construyendo ahora...${NC}"
    echo -e "${YELLOW}Este paso puede tardar varios minutos.${NC}"
    
    if docker build \
        --build-arg FASTA_URL_ARG="$FASTA_URL" \
        --build-arg GFF_URL_ARG="$GFF_URL" \
        -t "$DOCKER_IMAGE_TAG" .; then
        echo -e "   ${GREEN}[OK] Imagen Docker construida con éxito.${NC}"
    else
        echo -e "${RED}[ERROR] Falló la construcción de la imagen. Revisa los mensajes de error.${NC}"
        echo -e "${YELLOW}Intentando limpiar las capas de imagen fallidas...${NC}"
        docker image prune -f
        echo -e "${YELLOW}Limpieza completada. Abortando la instalación.${NC}"
        exit 1
    fi
else
    echo "   [OK] La imagen Docker '${DOCKER_IMAGE_TAG}' ya existe."
fi

# --- PASO 4: Crear Estructura de Directorios ---
echo -e "\n${GREEN}4. Creando la estructura de directorios...${NC}"
mkdir -p config data/raw_fastq results/{qc,trimmed,aligned,variants,annotated,logs,comparison}
echo "   [OK] Estructura de directorios creada."

# --- PASO 5: Configurar Muestras ---
echo -e "\n${GREEN}5. Configurando las muestras...${NC}"
declare -a SAMPLES_ARRAY
echo -e "${YELLOW}Ingresa los nombres de tus muestras (prefijos de los archivos FASTQ). Presiona Enter para terminar.${NC}"
echo -e "${YELLOW}IMPORTANTE: Ingresa primero la muestra 'control' (ej. susceptible) y luego la 'caso' (ej. resistente).${NC}"
while true; do read -p "Nombre de muestra: " sample_name; [ -z "$sample_name" ] && break; SAMPLES_ARRAY+=("$sample_name"); done
if [ ${#SAMPLES_ARRAY[@]} -lt 2 ]; then echo -e "${RED}Se necesitan al menos dos muestras para la comparación. Abortando.${NC}"; exit 1; fi
echo "sample" > config/samples.tsv; for s in "${SAMPLES_ARRAY[@]}"; do echo "$s" >> config/samples.tsv; done
cat > config/config.yaml <<- EOM
samples: "config/samples.tsv"
EOM
echo "   [OK] Archivos de configuración generados."

# --- PASO 6: Generar Scripts de Ejecución ---
echo -e "\n${GREEN}6. Creando scripts de ejecución personalizados...${NC}"
CORES_TOTAL=$(nproc 2>/dev/null || echo 8)
RAM_SUGGESTED=32 # Sugerir 32GB es un valor seguro y eficiente para GATK en este genoma

read -p "Número total de cores a usar [Sugerido: $CORES_TOTAL]: " USER_CORES; USER_CORES=${USER_CORES:-$CORES_TOTAL}
read -p "Memoria RAM para GATK HaplotypeCaller (GB) [Sugerido: ${RAM_SUGGESTED}]: " USER_RAM; USER_RAM=${USER_RAM:-$RAM_SUGGESTED}

cat > run_pipeline.sh <<- EOM
#!/bin/bash
echo "Iniciando pipeline LVAR con ${USER_CORES} cores y GATK HaplotypeCaller con ${USER_RAM}GB RAM..." >&2
docker run --rm -it \\
    --user "\$(id -u):\$(id -g)" \\
    -v "\$(pwd)/config:/pipeline/config" \\
    -v "\$(pwd)/data:/pipeline/data" \\
    -v "\$(pwd)/results:/pipeline/results" \\
    -v "\$(pwd)/Snakefile:/pipeline/Snakefile" \\
    "${DOCKER_IMAGE_TAG}" \\
    --cores ${USER_CORES} \\
    --config gatk_ram_gb=${USER_RAM} \\
    "\$@"
EOM
chmod +x run_pipeline.sh
echo "   [OK] Script de ejecución del pipeline './run_pipeline.sh' creado."

cat > run_in_container.sh <<- EOM
#!/bin/bash
COMMAND_TO_RUN="\${@:-bash}"
echo "Ejecutando en contenedor: '\${COMMAND_TO_RUN}'" >&2
docker run --rm -it \\
    --user "\$(id -u):\$(id -g)" \\
    -v "\$(pwd):/pipeline" \\
    -w "/pipeline" \\
    --entrypoint "" \\
    "${DOCKER_IMAGE_TAG}" \\
    \${COMMAND_TO_RUN}
EOM
chmod +x run_in_container.sh
echo "   [OK] Script de utilidad './run_in_container.sh' creado."

# --- PASO 7: Organizar Archivos FASTQ ---
echo -e "\n${GREEN}7. Organizando los archivos de datos FASTQ...${NC}"
read -e -p "Proporciona la RUTA ABSOLUTA a la carpeta con tus archivos .fastq.gz: " SOURCE_DIR
if [ -d "$SOURCE_DIR" ]; then
    for s in "${SAMPLES_ARRAY[@]}"; do
        cp "$SOURCE_DIR/${s}_R1.fastq.gz" "data/raw_fastq/" 2>/dev/null || echo -e "${RED}-> No se pudo copiar ${s}_R1.fastq.gz${NC}"
        cp "$SOURCE_DIR/${s}_R2.fastq.gz" "data/raw_fastq/" 2>/dev/null || echo -e "${RED}-> No se pudo copiar ${s}_R2.fastq.gz${NC}"
    done
else
    echo -e "${RED}Directorio no encontrado. Deberás mover los archivos FASTQ manualmente.${NC}"
fi

# --- Conclusión ---
echo -e "\n${GREEN}¡CONFIGURACIÓN COMPLETADA!${NC}"
echo "Para ejecutar el análisis completo y la comparación, usa:"
echo -e "${GREEN}./run_pipeline.sh${NC}"
echo "Para ejecutar comandos auxiliares (como bcftools) o explorar, usa:"
echo -e "${GREEN}./run_in_container.sh <tu_comando>${NC}"
