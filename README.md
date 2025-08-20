# LVAR: Análisis de Variantes en *Leishmania*
[![Snakemake](https://img.shields.io/badge/snakemake-core-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Docker](https://img.shields.io/badge/docker-engine-blue.svg)](https://www.docker.com/)

**LVAR** es un pipeline bioinformático optimizado para la identificación de variantes genéticas (SNPs e Indels) a partir de datos de secuenciación de genoma completo (WGS) de *Leishmania braziliensis*. Está diseñado para comparar aislados con diferentes fenotipos, como la susceptibilidad y resistencia al tratamiento leishmanicida.

Este proyecto está construido con un enfoque de **máxima reproducibilidad y facilidad de uso**. Todo el software, dependencias, genoma de referencia y bases de datos de anotación están encapsulados en una **única imagen de Docker**, eliminando la necesidad de instalar cualquier herramienta bioinformática en el sistema del usuario.

---

## Arquitectura y Reproducibilidad

El único requisito para el usuario es tener **Git** y **Docker**.

El `Dockerfile` del proyecto se encarga de:
1.  Instalar un entorno Conda basado en Mambaforge para evitar problemas de compatibilidad.
2.  Instalar todas las herramientas necesarias (`Snakemake`, `BWA`, `GATK4`, `SnpEff`, etc.).
3.  Descargar los archivos de referencia (genoma FASTA y anotación GFF) directamente desde TriTrypDB.
4.  Construir una base de datos de SnpEff.

Esto garantiza que cualquier usuario obtendrá exactamente los mismos resultados a partir de los mismos datos de entrada.

## Flujo de Trabajo del Pipeline

1.  **Control de Calidad (QC)**: `FastQC` y `MultiQC` para evaluar las lecturas crudas.
2.  **Limpieza de Lecturas (Trimming)**: `fastp` para eliminar adaptadores y bases de baja calidad.
3.  **Alineamiento**: `BWA-MEM` contra un genoma de referencia de *L. braziliensis* (ej. *L. braziliensis* MHOM/BR/75/M2904).
4.  **Post-procesamiento de BAM**: `GATK MarkDuplicates` y `Samtools` para ordenar e indexar.
5.  **Llamada de Variantes**: `GATK HaplotypeCaller` para identificar SNPs e Indels.
6.  **Anotación de Variantes**: `SnpEff` para predecir el impacto funcional de las variantes.

---

## Requisitos Previos

Solo necesitas tener instalados dos programas en tu sistema Linux:

-   **Git**: Para clonar el repositorio.
-   **Docker**: Para construir y ejecutar el entorno del pipeline.
    -   *Asegúrate de que el servicio de Docker esté activo (`sudo systemctl start docker`) y que tu usuario tenga permisos para ejecutarlo (o usa `sudo`).*

---

### Guia Completa

### **Fase 1: Configuración Inicial del Proyecto (Se realiza una sola vez)**

Esta fase prepara tu entorno de trabajo.

#### **1.1. Requisitos Previos**
Asegúrate de tener instalados los siguientes programas en tu sistema:
*   **Git**: Para clonar el repositorio.
*   **Docker**: Para construir y ejecutar el entorno del pipeline.

#### **1.2. Obtener el Código del Pipeline**
Clona el repositorio desde GitHub y entra en el directorio del proyecto.
```bash
git clone https://github.com/juanjoaguirre-IPE/lvar.git
cd lvar
```

#### **1.3. Ejecutar el Script de Configuración Interactivo**
El script `setup.sh` es un asistente que configurará todo por ti. Dale permisos de ejecución y lánzalo.
```bash
chmod +x setup.sh
./setup.sh
```
o:
```bash
bash setup.sh
```

El script te guiará a través de los siguientes pasos:

1.  **Configuración de Enlaces de Referencia**:
    *   El script te pedirá que introduzcas la URL para el genoma (`.fasta`) y la anotación (`.gff`).
    *   **Te ofrecerá un enlace por defecto**. Si el enlace sigue siendo válido, simplemente puedes **presionar Enter** para usarlo.
    *   Si el enlace está roto, puedes pegar uno nuevo. El script **verificará que el nuevo enlace sea válido** antes de continuar.

2.  **Construcción de la Imagen Docker**: Usando los enlaces que proporcionaste, construirá la imagen de Docker. Este proceso puede tardar más de 30 minutos la primera vez.

3.  **Configuración de Muestras**: Te pedirá que introduzcas los nombres de tus aislados (p. ej., `susceptible`, `resistente`).

4.  **Configuración de Recursos**: Te preguntará cuántos cores de CPU y memoria RAM (en GB) quieres asignar.

5.  **Organización de Datos**: Te pedirá la ruta a la carpeta donde tienes tus archivos `.fastq.gz` para copiarlos.

**Resultado de la Fase 1:** Tendrás una carpeta `lvar` con una estructura de proyecto completa y dos scripts clave: `run_pipeline.sh` (para el análisis principal) y `run_in_container.sh` (para tareas auxiliares).

---

### **Fase 2: Ejecución del Pipeline Principal**

Esta fase realiza el análisis bioinformático completo.

#### **2.1. Iniciar el Análisis**
Ejecuta el script `run_pipeline.sh`.
```bash
# Estando en la carpeta lvar
./run_pipeline.sh
```
Snakemake se iniciará, mostrará el plan de trabajo y comenzará el análisis.

#### **2.2. (Opcional) Realizar una Simulación (Dry-run)**
Para ver qué pasos se ejecutarían sin iniciar el análisis, usa el flag `-n`.
```bash
./run_pipeline.sh -n
```

**Resultado de la Fase 2:** La carpeta `results/` estará completamente poblada. Los archivos más importantes para el siguiente paso estarán en `results/annotated/`.

---

### **Fase 3: Análisis Comparativo Post-Pipeline**

Una vez que el pipeline principal ha terminado, encuentra las variantes únicas del aislado resistente.

#### **3.1. Preparar los Archivos VCF**
Usa el script auxiliar `run_in_container.sh` para ejecutar `bgzip` y `tabix`.
```bash
# Comprimir los VCFs
./run_in_container.sh bgzip results/annotated/susceptible.ann.vcf
./run_in_container.sh bgzip results/annotated/resistant.ann.vcf

# Indexar los VCFs comprimidos
./run_in_container.sh tabix -p vcf results/annotated/susceptible.ann.vcf.gz
./run_in_container.sh tabix -p vcf results/annotated/resistant.ann.vcf.gz
```

#### **3.2. Comparar los VCFs para Encontrar Variantes Únicas**
Usa `bcftools isec` para crear archivos con las diferencias.
```bash
./run_in_container.sh bcftools isec -p results/annotated/comparison \
    results/annotated/susceptible.ann.vcf.gz \
    results/annotated/resistant.ann.vcf.gz
```

**Resultado de la Fase 3:** Se creará una nueva carpeta `results/annotated/comparison/`. El archivo clave es **`0001.vcf`**, que contiene las variantes **exclusivas del aislado resistente**.

---

### **Fase 4: Visualización e Interpretación con IGV (en Windows/Mac/Linux)**

Este paso te permite ver tus hallazgos en su contexto genómico.

#### **4.1. Instalar IGV**
Descarga e instala el **Integrative Genomics Viewer (IGV)** en tu máquina local desde [aquí](https://software.broadinstitute.org/software/igv/download).

#### **4.2. Extraer los Archivos de Referencia del Contenedor**
Para que IGV entienda tus datos, necesitas el genoma y la anotación. Ejecuta estos comandos en tu carpeta `lvar` para extraerlos del contenedor.
```bash
# Extraer el genoma FASTA
./run_in_container.sh cat /opt/db/genome.fasta > leishmania_genome.fasta

# Extraer la anotación GFF
./run_in_container.sh cat /opt/conda/share/snpeff-5.2-1/data/Lbraziliensis_2019_manual/genes.gff > leishmania_annotation.gff
```
Ahora tendrás dos nuevos archivos en tu carpeta `lvar`.

#### **4.3. Cargar Datos en IGV**
Abre IGV y carga los archivos en este orden:
1.  **Cargar Genoma**: `Genomes` -> `Load Genome from File...` -> selecciona `leishmania_genome.fasta`.
2.  **Cargar Anotación**: `File` -> `Load from File...` -> selecciona `leishmania_annotation.gff`.
3.  **Cargar Variantes Candidatas**: `File` -> `Load from File...` -> selecciona `results/annotated/comparison/0001.vcf`.

#### **4.4. Investigar Genes Candidatos**
Puedes obtener las coordenadas de los genes de interés en Tritrypdb. Utiliza esas coordenadas en la barra de búsqueda en la parte superior de IGV. Observa la pista de tu archivo `0001.vcf`. Si ves una línea de color sobre la estructura de un gen, has encontrado una mutación candidata. Haz clic en ella para ver los detalles de la anotación de SnpEff.
