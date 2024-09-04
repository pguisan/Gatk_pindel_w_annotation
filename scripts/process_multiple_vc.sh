#!/bin/bash

# Ruta completa a mamba
MAMBA_PATH="/home/paula/conda/bin/mamba"

# Ruta al entorno gatk_env
GATK_ENV="/home/paula/conda/envs/gatk_env"

# Verificar si el entorno gatk_env existe
if [ ! -d "$GATK_ENV" ]; then
    echo "Error: El entorno gatk_env no existe en $GATK_ENV. Verifica la ruta."
    exit 1
fi

echo "Usando entorno gatk_env en: $GATK_ENV"

# Definir las rutas
SCRIPTS_DIR="$(pwd)"
BASE_DIR="$(dirname "$SCRIPTS_DIR")"
SAMPLE_GROUPS=("IMR" "CZA")

echo "Directorio base: $BASE_DIR"
echo "Directorio de scripts: $SCRIPTS_DIR"

# Función para procesar cada muestra
process_sample() {
    local SAMPLE_DIR=$1
    local SAMPLE_NAME=$(basename "$SAMPLE_DIR")
    echo "Procesando muestra: $SAMPLE_NAME en $SAMPLE_DIR"

    # Encontrar el directorio de lecturas
    local READS_DIR=$(find "$SAMPLE_DIR" -type d -name "*_reads" | head -n 1)
    if [ -z "$READS_DIR" ]; then
        echo "Error: No se encontró directorio de lecturas para $SAMPLE_NAME"
        return 1
    fi

    # Encontrar el directorio indexado y los archivos de referencia
    local INDEXADO_DIR=$(find "$SAMPLE_DIR" -type d -name "*_indexado" | head -n 1)
    local REF_FASTA=$(find "$INDEXADO_DIR" -type f -name "*.fasta" | head -n 1)
    local REF_GBFF=$(find "$INDEXADO_DIR" -type f -name "*.gbff" | head -n 1)

    if [ -z "$REF_FASTA" ] || [ -z "$REF_GBFF" ]; then
        echo "Error: No se encontraron archivos de referencia para $SAMPLE_NAME"
        return 1
    fi

    # Crear el índice BWA si no existe
    if [ ! -f "${REF_FASTA}.bwt" ]; then
        echo "Creando índice BWA..."
        bwa index $REF_FASTA
    fi

    # Crear el índice del genoma de referencia si no existe
    if [ ! -f "${REF_FASTA}.fai" ]; then
        echo "Creando índice del genoma de referencia..."
        samtools faidx $REF_FASTA
    fi

    # Crear el diccionario de secuencias si no existe
    if [ ! -f "${REF_FASTA%.fasta}.dict" ]; then
        echo "Creando diccionario de secuencias..."
        picard CreateSequenceDictionary R=$REF_FASTA O=${REF_FASTA%.fasta}.dict
    fi

    # Procesar cada subdirectorio en el directorio de lecturas
    for SUB_SAMPLE_DIR in "$READS_DIR"/*; do
        if [ -d "$SUB_SAMPLE_DIR" ]; then
            local SUB_SAMPLE_NAME=$(basename "$SUB_SAMPLE_DIR")
            echo "Procesando sub-muestra: $SUB_SAMPLE_NAME"

            # Crear directorio para los resultados de la sub-muestra
            local RESULTS_DIR="$SAMPLE_DIR/${SUB_SAMPLE_NAME}_results"
            mkdir -p "$RESULTS_DIR"
            echo "Directorio de resultados creado: $RESULTS_DIR"

            # Encontrar y renombrar los archivos de lectura si es necesario
            local R1_FILE=$(find "$SUB_SAMPLE_DIR" -name "*_R1_001.fastq.gz" | head -n 1)
            local R2_FILE=$(find "$SUB_SAMPLE_DIR" -name "*_R2_001.fastq.gz" | head -n 1)

            if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
                for READ_FILE in "$SUB_SAMPLE_DIR"/*.f*q.gz; do
                    if [[ $READ_FILE =~ .*_R1.*\.f.*q\.gz$ ]]; then
                        NEW_NAME="$SUB_SAMPLE_DIR/${SUB_SAMPLE_NAME}_R1_001.fastq.gz"
                        mv "$READ_FILE" "$NEW_NAME"
                        R1_FILE="$NEW_NAME"
                        echo "Renombrado: $READ_FILE -> $R1_FILE"
                    elif [[ $READ_FILE =~ .*_R2.*\.f.*q\.gz$ ]]; then
                        NEW_NAME="$SUB_SAMPLE_DIR/${SUB_SAMPLE_NAME}_R2_001.fastq.gz"
                        mv "$READ_FILE" "$NEW_NAME"
                        R2_FILE="$NEW_NAME"
                        echo "Renombrado: $READ_FILE -> $R2_FILE"
                    fi
                done
            fi

            if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
                echo "Error: No se encontraron archivos de lectura para $SUB_SAMPLE_NAME"
                echo "Contenido del directorio de lecturas:"
                ls -l "$SUB_SAMPLE_DIR"
                continue
            fi

            echo "Archivos encontrados para $SUB_SAMPLE_NAME:"
            echo "Referencia FASTA: $REF_FASTA"
            echo "Referencia GBFF: $REF_GBFF"
            echo "Lecturas R1: $R1_FILE"
            echo "Lecturas R2: $R2_FILE"

            # Alineamiento
            echo "Realizando alineamiento para $SUB_SAMPLE_NAME"
            bwa mem $REF_FASTA $R1_FILE $R2_FILE > "$RESULTS_DIR/${SUB_SAMPLE_NAME}_aligned.sam"

            # Convertir SAM a BAM
            echo "Convirtiendo SAM a BAM para $SUB_SAMPLE_NAME"
            samtools view -bS "$RESULTS_DIR/${SUB_SAMPLE_NAME}_aligned.sam" > "$RESULTS_DIR/${SUB_SAMPLE_NAME}_aligned.bam"

            # Ordenar BAM
            echo "Ordenando BAM para $SUB_SAMPLE_NAME"
            samtools sort "$RESULTS_DIR/${SUB_SAMPLE_NAME}_aligned.bam" -o "$RESULTS_DIR/${SUB_SAMPLE_NAME}_sorted.bam"

            # Indexar BAM ordenado
            echo "Indexando BAM ordenado para $SUB_SAMPLE_NAME"
            samtools index "$RESULTS_DIR/${SUB_SAMPLE_NAME}_sorted.bam"

            # Añadir grupos de lectura
            echo "Añadiendo grupos de lectura para $SUB_SAMPLE_NAME"
            picard AddOrReplaceReadGroups \
                I="$RESULTS_DIR/${SUB_SAMPLE_NAME}_sorted.bam" \
                O="$RESULTS_DIR/${SUB_SAMPLE_NAME}_with_RG.bam" \
                RGID=4 \
                RGLB=lib1 \
                RGPL=ILLUMINA \
                RGPU=unit1 \
                RGSM=$SUB_SAMPLE_NAME \
                VALIDATION_STRINGENCY=LENIENT

            # Marcar duplicados
            echo "Marcando duplicados para $SUB_SAMPLE_NAME"
            picard MarkDuplicates \
                I="$RESULTS_DIR/${SUB_SAMPLE_NAME}_with_RG.bam" \
                O="$RESULTS_DIR/${SUB_SAMPLE_NAME}_deduped.bam" \
                M="$RESULTS_DIR/${SUB_SAMPLE_NAME}_metrics.txt" \
                VALIDATION_STRINGENCY=LENIENT

            # Indexar BAM final
            echo "Indexando BAM final para $SUB_SAMPLE_NAME"
            samtools index "$RESULTS_DIR/${SUB_SAMPLE_NAME}_deduped.bam"

            # HaplotypeCaller (strict)
            echo "Ejecutando HaplotypeCaller para $SUB_SAMPLE_NAME"
            gatk HaplotypeCaller \
                -R $REF_FASTA \
                -I "$RESULTS_DIR/${SUB_SAMPLE_NAME}_deduped.bam" \
                -O "$RESULTS_DIR/${SUB_SAMPLE_NAME}_output_strict.vcf.gz" \
                --sample-ploidy 1 \
                --min-base-quality-score 20 \
                --standard-min-confidence-threshold-for-calling 30 \
                --min-pruning 2 \
                --min-dangling-branch-length 2 \
                --annotation-group StandardAnnotation \
                --annotation ReadPosRankSumTest \
                --annotation QualByDepth \
                --annotation FisherStrand

            # VariantFiltration
            echo "Ejecutando VariantFiltration para $SUB_SAMPLE_NAME"
            gatk VariantFiltration \
                -R $REF_FASTA \
                -V "$RESULTS_DIR/${SUB_SAMPLE_NAME}_output_strict.vcf.gz" \
                -O "$RESULTS_DIR/${SUB_SAMPLE_NAME}_filtered_strict.vcf.gz" \
                --filter-expression "DP < 4 || QD < 2.0 || FS > 60.0 || MQ < 30.0 || SOR > 4.0" \
                --filter-name "BalancedFilter"

            # Ejecutar el script de Python para anotación de variantes de GATK
            echo "Anotando resultados de GATK para la sub-muestra $SUB_SAMPLE_NAME"
            python $SCRIPTS_DIR/vcf_annotator_gatk.py \
                "$RESULTS_DIR/${SUB_SAMPLE_NAME}_filtered_strict.vcf.gz" \
                $REF_GBFF \
                "$RESULTS_DIR/${SUB_SAMPLE_NAME}_gatk_annotated.tsv"

            # Ejecutar Pindel para indels
            echo "Ejecutando Pindel para la sub-muestra $SUB_SAMPLE_NAME"
            # Crear el archivo de configuración de Pindel
            echo -e "$RESULTS_DIR/${SUB_SAMPLE_NAME}_sorted.bam\t400\t$SUB_SAMPLE_NAME" > "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_config.txt"
            
            pindel \
                -f $REF_FASTA \
                -i "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_config.txt" \
                -c ALL \
                -o "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel" \
                -x 5 \
                -M 3 \
                -A 30 \
                -T 12

            # Convertir salida de Pindel a VCF
            echo "Convirtiendo salida de Pindel a VCF para $SUB_SAMPLE_NAME"
            pindel2vcf \
                -p "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_D" \
                -r $REF_FASTA \
                -R "Klebsiella_pneumoniae_reference" \
                -d 20200101 \
                -v "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_indels.vcf" \
				-G

            # Anotar los resultados de indels de Pindel
            echo "Anotando resultados de indels de Pindel para la sub-muestra $SUB_SAMPLE_NAME"
            python $SCRIPTS_DIR/vcf_annotator_pindel.py \
                "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_indels.vcf" \
                "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_indel_annotated.tsv" \
                $REF_GBFF

            echo "Anotación de indels de Pindel completada para la sub-muestra $SUB_SAMPLE_NAME"

            echo "Procesamiento completado para la sub-muestra $SUB_SAMPLE_NAME"
            
            # Llamar a la función de limpieza después de procesar la sub-muestra
            cleanup_intermediate_files "$RESULTS_DIR" "$SUB_SAMPLE_NAME"
        fi
    done
}
cleanup_intermediate_files() {
    local RESULTS_DIR=$1
    local SUB_SAMPLE_NAME=$2

    echo "Limpiando archivos intermedios para $SUB_SAMPLE_NAME en $RESULTS_DIR"

    # Eliminar archivos SAM y BAM intermedios
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_aligned.sam"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_aligned.bam"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_sorted.bam"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_sorted.bam.bai"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_with_RG.bam"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_deduped.bam"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_deduped.bam.bai"

    # Eliminar archivos VCF intermedios
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_output_strict.vcf.gz"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_output_strict.vcf.gz.tbi"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_filtered_strict.vcf.gz.tbi"

    # Eliminar archivos de métricas y otros archivos intermedios
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_metrics.txt"

    # Eliminar archivos intermedios de Pindel
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_config.txt"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_D"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_SI"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_LI"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_INV"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_TD"
    rm -f "$RESULTS_DIR/${SUB_SAMPLE_NAME}_pindel_RP"
    # Eliminar archivos de Pindel de 0KB
    find "$RESULTS_DIR" -name "${SUB_SAMPLE_NAME}_pindel_*" -size 0 -delete

    echo "Limpieza completada para $SUB_SAMPLE_NAME"

    # Verificar si se eliminaron los archivos
    remaining_files=$(find "$RESULTS_DIR" -type f | grep -v "_gatk_annotated.tsv" | grep -v "_pindel_indel_annotated.tsv")
    if [ -n "$remaining_files" ]; then
        echo "Advertencia: Algunos archivos intermedios no se eliminaron:"
        echo "$remaining_files"
    else
        echo "Todos los archivos intermedios fueron eliminados correctamente."
    fi
}

# Recorrer todas las muestras en IMR y CZA
for GROUP in "${SAMPLE_GROUPS[@]}"; do
    GROUP_DIR="$BASE_DIR/$GROUP"
    echo "Buscando muestras en: $GROUP_DIR"
    
    if [ ! -d "$GROUP_DIR" ]; then
        echo "Advertencia: El directorio $GROUP_DIR no existe."
        continue
    fi
    
    SAMPLE_COUNT=0
    for SAMPLE_DIR in "$GROUP_DIR"/*; do
        if [ -d "$SAMPLE_DIR" ]; then
            echo "Encontrada muestra: $(basename "$SAMPLE_DIR")"
            process_sample "$SAMPLE_DIR"
            SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
        fi
    done
    
    echo "Procesadas $SAMPLE_COUNT muestras en $GROUP"
done

echo "Procesamiento completo para todas las muestras."

