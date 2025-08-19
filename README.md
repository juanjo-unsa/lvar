# LVAR: Análisis de Variantes en Leishmania
[![Snakemake](https://img.shields.io/badge/snakemake-core-brightgreen.svg)](https://snakemake.readthedocs.io)
[![Docker](https://img.shields.io/badge/docker-engine-blue.svg)](https://www.docker.com/)

**LVAR** es un pipeline bioinformático robusto y reproducible para la identificación de variantes genéticas (SNPs e Indels) a partir de datos de secuenciación de genoma completo (WGS) de *Leishmania braziliensis*. Está diseñado para comparar aislados con diferentes fenotipos, como la susceptibilidad y resistencia a tratamientos farmacológicos.

Este proyecto está construido con un enfoque de **máxima reproducibilidad y facilidad de uso**. Todo el software, dependencias, genoma de referencia y bases de datos de anotación están encapsulados en una **única imagen de Docker**, eliminando la necesidad de instalar cualquier herramienta bioinformática en el sistema del usuario.

---

## Arquitectura y Reproducibilidad

El único requisito para el usuario es tener **Git** y **Docker**.

El `Dockerfile` del proyecto se encarga de:
1.  Instalar un entorno Conda basado en Mambaforge para evitar problemas de compatibilidad.
2.  Instalar todas las herramientas necesarias (`Snakemake`, `BWA`, `GATK4`, `SnpEff`, etc.).
3.  Descargar los archivos de referencia (genoma FASTA y anotación GFF) directamente desde TriTrypDB.
4.  Construir una base de datos de SnpEff.

Esto garantiza que cualquier persona, en cualquier máquina con Docker, obtendrá exactamente los mismos resultados a partir de los mismos datos de entrada.

## Flujo de Trabajo del Pipeline

1.  **Control de Calidad (QC)**: `FastQC` y `MultiQC` para evaluar las lecturas crudas.
2.  **Limpieza de Lecturas (Trimming)**: `fastp` para eliminar adaptadores y bases de baja calidad.
3.  **Alineamiento**: `BWA-MEM` contra el genoma de referencia de *L. braziliensis* (MHOM/BR/75/M2904).
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

## Guía de Inicio Rápido

Sigue estos tres sencillos pasos para poner en marcha el proyecto.

### 1. Clonar el Repositorio

Abre una terminal y clona este repositorio en tu máquina. Luego, entra en el directorio recién creado.
```bash
git clone https://github.com/juanjo-unsa/lvar.git
cd lvar
```
### 2. Instalación
```
chmod +x setup.sh
./setup.sh
```
### 3. Ejecutar el pipeline
```
./run_pipeline.sh
```
