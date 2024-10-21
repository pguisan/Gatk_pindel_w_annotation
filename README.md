# 🧬 GATK-Pindel Genomic Analysis Pipeline

## 🎓 Doctoral Thesis Project

This project is part of a doctoral thesis in bioinformatics focused on the analysis of antibiotic resistance in pathogenic bacteria.

👩‍🔬 **PhD candidate**: Paula Guijarro-Sánchez

🏆 Xunta de Galicia Predoctoral Student Grant – IN606A- 2021/021

### 👥 Supervisors:

- 🩺 Dr. Alejandro Beceiro Casas (SERGAS)
- 💻 Dr. Carlos Fernandez-Lozano (UDC)

## 🏆 Featured Publication

The results obtained with this pipeline have been published in a prestigious scientific journal:

📚 **Journal**: International Journal of Antimicrobial Agents (Q1 in Microbiology)
🔗 **DOI**: 10.1016/j.ijantimicag.2023.106935

👩‍🔬 **Co-first author**: Paula Guijarro-Sánchez

This publication emphasizes the critical importance of selecting the appropriate combination of variant calling tools to effectively identify mutations that link antibiotic resistance phenotypes with genomic changes in *Klebsiella pneumoniae*. It demonstrates how the careful integration of different variant calling approaches can enhance the detection of resistance-associated mutations, providing valuable insights into the genomic basis of antibiotic resistance.

## 🎯 Purpose of the Analysis

This bioinformatics pipeline was developed for the analysis of genomic variants associated with antibiotic resistance in Klebsiella pneumoniae. Specifically, the study focuses on:

1. 💊 Analysis of resistance to imipenem/relebactam and ceftazidime/avibactam.
2. 🔬 Detection of genomic variants in *K. pneumoniae* strains exposed to increasing concentrations of these antibiotics.
3. 🧫 Identification of mutations associated with resistance, especially in genes such as *bla*KPC and outer membrane proteins.

The pipeline processes Illumina sequencing data from parental strains and selected resistant mutants, allowing the identification of genomic variants that could be associated with observed resistance phenotypes.

---

This repository contains a bioinformatics pipeline for the analysis of genomic variants in bacteria, with a particular focus on detecting antibiotic resistance mechanisms. The pipeline integrates cutting-edge tools such as GATK and Pindel for comprehensive variant analysis.

## ✨ Main Features

-🧹 Preprocessing of sequencing data

-🔍 Read alignment with BWA-MEM

-🧩 Detection of SNVs and small indels with GATK4 HaplotypeCaller

-📏 Detection of large indels with Pindel

-📝 Functional annotation of variants

-📊 Generation of detailed reports in TSV format

## 🛠 Requirements

- 🐚 Bash
- 🐍 Mamba/Conda
- 🧬 BWA
- 🧰 Samtools
- 🃏 Picard
- 🧬 GATK4
- 📏 Pindel
- 🐍 Python 3.x

## 📂 Repository Structure
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

Within the CZA and IMR folders, the following structure is expected for each sample:

```
📁[Nombre de la muestra]/
├── 📁[Nombre de la muestra]_reads/
│   ├── 📄[archivos de lecturas .fastq.gz]
└── 📁[Nombre de la muestra]_indexado/
    ├── 📄[archivos del genoma indexado]
```

## 🚀 Usage

1. Clone the repository:
   ```
   git clone https://github.com/pguisan/Gatk_pindel_w_annotation.git
   ```

2. Place your sequencing data and indexed genomes in the corresponding folders within `CZA` or `IMR`.

3. Run the main script:
   ```
   bash scripts/process_multiple_vc.sh [ruta al directorio de la muestra]
   ```

4. Results will be stored in specific folders for each sample within  `CZA` or `IMR`.

## 🔬 Analysis Pipeline

The pipeline follows these main steps:

1. 🧹 Data preprocessing
2. 🔍 Alignment with BWA-MEM
3. 🔧 Post-alignment processing (SAM to BAM conversion, sorting, indexing)
4. 🧬 Variant calling with GATK4 HaplotypeCaller
5. 🧫 Variant filtering with GATK4 VariantFiltration
6. 📏 Large indel detection with Pindel
7. 📝 Functional annotation of variants

## ⚙️ Customization

You can adjust the parameters of each tool by editing the corresponding scripts in the  `scripts/` folder. 

## ✅ Validation

The pipeline has been validated by comparing its results with those obtained using Snippy, an established tool for variant detection in bacterial genomes.

## 👥 Contributions

Contributions are welcome. Please open an issue to discuss major changes before submitting a pull request.
