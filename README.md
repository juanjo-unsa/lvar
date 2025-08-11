# LVAR: Análisis de Variantes en Leishmania
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Docker](https://img.shields.io/badge/docker-engine-blue.svg)](https://www.docker.com/)

**LVAR** es un pipeline bioinformático para la identificación de variantes genéticas (SNPs e Indels) a partir de datos de secuenciación de genoma completo (WGS) de *Leishmania braziliensis*. Está diseñado específicamente para comparar aislados con diferentes fenotipos.

El pipeline está construido con **Snakemake** para la gestión del flujo de trabajo y utiliza **Docker** para encapsular cada herramienta bioinformática, garantizando una ejecución idéntica en cualquier sistema Linux.

---

## Características Principales

-   **Reproducibilidad Total**: Gracias a Snakemake y Docker, los resultados son consistentes en cualquier máquina.
-   **Instalación Sencilla**: Un único script (`install.sh`) configura el entorno, verifica dependencias y prepara el proyecto.
-   **Configuración de Recursos**: El script de instalación te permite adaptar el uso de CPU (cores) y RAM a los recursos de tu servidor.
-   **Análisis Completo**: El pipeline abarca desde el control de calidad de las lecturas crudas hasta la anotación funcional de las variantes.

## Flujo de Trabajo del Pipeline (Workflow)

1.  **Control de Calidad (QC)**: Se evalúa la calidad de las lecturas FASTQ crudas usando `FastQC` y se genera un reporte agregado con `MultiQC`.
2.  **Limpieza de Lecturas (Trimming)**: Se eliminan adaptadores y bases de baja calidad con `fastp`.
3.  **Alineamiento**: Las lecturas limpias se alinean contra un genoma de referencia de *L. braziliensis* usando `BWA-MEM`.
4.  **Post-procesamiento de BAM**: Los alineamientos se ordenan, y los duplicados de PCR se marcan con `Samtools` y `GATK MarkDuplicates`.
5.  **Llamada de Variantes (Variant Calling)**: Se identifican SNPs e Indels para cada muestra individualmente usando `GATK HaplotypeCaller`.
6.  **Anotación de Variantes**: Se predice el impacto funcional de las variantes (ej. missense, frameshift, stop-gained) con `SnpEff`.
---

## Requisitos Previos

Antes de comenzar, asegúrate de tener instalados los siguientes programas en tu sistema Linux:

-   **Git**: Para clonar el repositorio.
-   **Docker**: Para ejecutar las herramientas bioinformáticas en contenedores.
    -   *Asegúrate de que el servicio de Docker esté en ejecución (`sudo systemctl start docker`).*
-   **Snakemake**: Para orquestar el pipeline. Se recomienda instalarlo a través de Conda:
    ```bash
    conda install -c bioconda snakemake
    ```
---

## Guía de Inicio Rápido

El único requisito de software en tu sistema es **Git** y **Docker**.

### 1. Clonar el Repositorio

Abre una terminal y clona este repositorio en tu máquina:
```bash
git clone https://github.com/juanjo-unsa/lvar.git
cd lvar
