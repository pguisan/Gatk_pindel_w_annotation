# 🧬 Pipeline de Análisis Genómico GATK-Pindel

## 🎓 Proyecto de Tesis Doctoral

Este proyecto es parte de una tesis doctoral en bioinformática centrada en el análisis de resistencia a antibióticos en bacterias patógenas.

👩‍🔬 **PhD candidate:** Paula Guijarro-Sánchez
🏆 Xunta de Galicia Predoctoral Student Grant – IN606A- 2021/021

### 👥 Supervisores:

- 🩺 Dr. Alejandro Beceiro Casas (SERGAS)
- 💻 Dr. Carlos Fernandez-Lozano (UDC)

## 🏆 Publicación Destacada

Los resultados obtenidos con este pipeline han sido publicados en una prestigiosa revista científica:

📚 **Revista:** International Journal of Antimicrobial Agents (Q1 en Microbiología)
🔗 **DOI:** [10.1016/j.ijantimicag.2023.106935](https://doi.org/10.1016/j.ijantimicag.2023.106935)
👩‍🔬 **Co-primera autora:** Paula Guijarro-Sánchez

Esta publicación destaca la importancia y el impacto del trabajo realizado, validando la eficacia del pipeline desarrollado en el análisis de resistencia a antibióticos en *Klebsiella pneumoniae*.


## 🎯 Propósito del Análisis

Este pipeline bioinformático se desarrolló para el análisis de variantes genómicas asociadas con la resistencia a antibióticos en *Klebsiella pneumoniae*. Específicamente, el estudio se centra en:

1. 💊 Análisis de resistencia a imipenem/relebactam y ceftazidima/avibactam.
2. 🔬 Detección de variantes genómicas en cepas de *K. pneumoniae* expuestas a concentraciones crecientes de estos antibióticos.
3. 🧫 Identificación de mutaciones asociadas con la resistencia, especialmente en genes como *bla*KPC y proteínas de la membrana externa.

El pipeline procesa datos de secuenciación de Illumina de cepas parentales y mutantes resistentes seleccionados, permitiendo la identificación de variantes genómicas que podrían estar asociadas con fenotipos de resistencia observados.

---

Este repositorio contiene un pipeline bioinformático para el análisis de variantes genómicas en bacterias, con un enfoque particular en la detección de mecanismos de resistencia a antibióticos. El pipeline integra herramientas de vanguardia como GATK y Pindel para un análisis exhaustivo de variantes.

## 📋 Tabla de Contenidos

- [Características Principales](#-características-principales)
- [Requisitos](#-requisitos)
- [Estructura del Repositorio](#-estructura-del-repositorio)
- [Uso](#-uso)
- [Pipeline de Análisis](#-pipeline-de-análisis)
- [Personalización](#-personalización)
- [Validación](#-validación)
- [Contribuciones](#-contribuciones)

## ✨ Características Principales

- 🧹 Preprocesamiento de datos de secuenciación
- 🔍 Alineamiento de lecturas con BWA-MEM
- 🧩 Detección de SNVs e indels pequeños con GATK4 HaplotypeCaller
- 📏 Detección de indels grandes con Pindel
- 📝 Anotación funcional de variantes
- 📊 Generación de informes detallados en formato TSV

## 🛠 Requisitos

- 🐚 Bash
- 🐍 Mamba/Conda
- 🧬 BWA
- 🧰 Samtools
- 🃏 Picard
- 🧬 GATK4
- 📏 Pindel
- 🐍 Python 3.x

## 📂 Estructura del Repositorio
```
./
├── 📁scripts/
│   ├── 📜process_multiple_vc.sh
│   ├── 📜vcf_annotator_gatk.py
│   └── 📜vcf_annotator_pindel.py
├── 📁CZA/
│   └── [carpeta vacía para muestras de CZA]
├── 📁IMR/
│   └── [carpeta vacía para muestras de IMR]
└── 📄README.md
```

Dentro de las carpetas `CZA` e `IMR`, se espera la siguiente estructura para cada muestra:

```
📁[Nombre de la muestra]/
├── 📁[Nombre de la muestra]_reads/
│   ├── 📄[archivos de lecturas .fastq.gz]
└── 📁[Nombre de la muestra]_indexado/
    ├── 📄[archivos del genoma indexado]
```

## 🚀 Uso

1. Clone el repositorio:
   ```
   git clone https://github.com/pguisan/Gatk_pindel_w_annotation.git
   ```

2. Coloque sus datos de secuenciación y genomas indexados en las carpetas correspondientes dentro de `CZA` o `IMR`.

3. Ejecute el script principal:
   ```
   bash scripts/process_multiple_vc.sh [ruta al directorio de la muestra]
   ```

4. Los resultados se almacenarán en carpetas específicas para cada muestra dentro de `CZA` o `IMR`.
4. Los resultados se almacenarán en carpetas específicas para cada muestra dentro de `CZA` o `IMR`.

## 🔬 Pipeline de Análisis

El pipeline sigue estos pasos principales:

1. 🧹 Preprocesamiento de datos
2. 🔍 Alineamiento con BWA-MEM
3. 🔧 Procesamiento post-alineamiento (conversión SAM a BAM, ordenamiento, indexación)
4. 🧬 Llamado de variantes con GATK4 HaplotypeCaller
5. 🧫 Filtrado de variantes con GATK4 VariantFiltration
6. 📏 Detección de indels grandes con Pindel
7. 📝 Anotación funcional de variantes

## ⚙️ Personalización

Puede ajustar los parámetros de cada herramienta editando los scripts correspondientes en la carpeta `scripts/`. 

## ✅ Validación

El pipeline ha sido validado comparando sus resultados con los obtenidos usando Snippy, una herramienta establecida para la detección de variantes en genomas bacterianos.

## 👥 Contribuciones

Las contribuciones son bienvenidas. Por favor, abra un issue para discutir cambios mayores antes de enviar un pull request.
