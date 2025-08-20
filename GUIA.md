### Protocolo Completo: Análisis de Variantes en *Leishmania* con el Pipeline LVAR

#### **Resumen del Protocolo**
Este documento describe el procedimiento completo para configurar y ejecutar el pipeline LVAR, desde la obtención del código hasta el análisis comparativo de los resultados. El pipeline está diseñado para ser ejecutado en un entorno Linux con Docker y Git instalados.

---

### **Fase 1: Configuración Inicial del Proyecto (Se realiza una sola vez)**

Esta fase prepara tu entorno de trabajo.

#### **1.1. Requisitos Previos**
Asegúrate de tener instalados los siguientes programas en tu sistema:
*   **Git**: Para clonar el repositorio.
    *   Verificar instalación: `git --version`
*   **Docker**: Para construir y ejecutar el entorno del pipeline.
    *   Verificar instalación: `docker --version`
    *   **Importante**: Asegúrate de que el servicio de Docker esté activo (`sudo systemctl start docker`) y que tu usuario tenga permisos para ejecutarlo.

#### **1.2. Obtener el Código del Pipeline**
Clona el repositorio desde GitHub y entra en el directorio del proyecto.
```bash
git clone https://github.com/juanjo-unsa/lvar.git
cd lvar
```

#### **1.3. Ejecutar el Script de Configuración Interactivo**
El script `setup.sh` es un asistente que configurará todo por ti. Dale permisos de ejecución y lánzalo.
```bash
chmod +x setup.sh
./setup.sh
```

El script te guiará a través de los siguientes pasos:
1.  **Construcción de la Imagen Docker**: Esto puede tardar más de 30 minutos la primera vez. Es el paso más largo, ya que descarga todas las herramientas y los datos de referencia.
2.  **Configuración de Muestras**: Te pedirá que introduzcas los nombres de tus aislados (p. ej., `susceptible`, `resistente`). Estos deben coincidir con el prefijo de tus archivos FASTQ.
3.  **Configuración de Recursos**: Te preguntará cuántos cores de CPU y cuánta memoria RAM (en GB) quieres asignar al pipeline.
4.  **Organización de Datos**: Finalmente, te pedirá la ruta a la carpeta donde tienes guardados tus archivos `.fastq.gz` y los copiará al lugar correcto.

**Resultado de la Fase 1:** Tendrás una carpeta `lvar` que contiene una estructura de proyecto completa, con tus datos en `data/raw_fastq/` y dos scripts clave: `run_pipeline.sh` y `run_in_container.sh`.

---

### **Fase 2: Ejecución del Pipeline Principal**

Esta fase realiza el análisis bioinformático desde los archivos FASTQ crudos hasta los VCF anotados.

#### **2.1. Iniciar el Análisis**
Ejecuta el script `run_pipeline.sh` que se creó en la fase anterior.
```bash
# Estando en la carpeta lvar
./run_pipeline.sh
```
Snakemake se iniciará, mostrará el plan de los trabajos a ejecutar y comenzará el análisis. Este proceso puede tardar varias horas o días dependiendo del tamaño de tus datos y los recursos del servidor.

#### **2.2. (Opcional) Realizar una Simulación (Dry-run)**
Para ver qué pasos se ejecutarían sin iniciar el análisis, usa el flag `-n`.
```bash
./run_pipeline.sh -n
```

**Resultado de la Fase 2:** La carpeta `results/` estará completamente poblada. Los archivos más importantes para el siguiente paso estarán en `results/annotated/`.

---

### **Fase 3: Análisis Comparativo Post-Pipeline**

Una vez que el pipeline principal ha terminado, el siguiente paso es encontrar las variantes únicas del aislado resistente.

#### **3.1. Preparar los Archivos VCF**
Usaremos el script auxiliar `run_in_container.sh` para ejecutar `bgzip` y `tabix`.
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

**Resultado de la Fase 3:** Se creará una nueva carpeta `results/annotated/comparison/`. El archivo clave dentro de esta carpeta es **`0001.vcf`**, que contiene las variantes que son **exclusivas del aislado resistente**.

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
1.  **Cargar Genoma**: Ve a `Genomes` -> `Load Genome from File...` y selecciona `leishmania_genome.fasta`.
2.  **Cargar Anotación**: Ve a `File` -> `Load from File...` y selecciona `leishmania_annotation.gff`. Verás las estructuras de los genes como cajas azules.
3.  **Cargar Variantes Candidatas**: Ve a `File` -> `Load from File...` y selecciona `results/annotated/comparison/0001.vcf`.

#### **4.4. Investigar Genes Candidatos**
Busca en Tritrypdb los genes de interes y obtiene sus coordenadas.
En IGV ubícate en las coordenadas y observa la pista de tu archivo `0001.vcf`. Si ves una línea de color sobre la estructura de un gen, has encontrado una mutación candidata. Haz clic en ella para ver los detalles de la anotación de SnpEff en el panel inferior.
