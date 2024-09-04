#!/usr/bin/env python3
import sys
import os
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from cyvcf2 import VCF
import csv
from datetime import datetime

# Configurar logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Obtener la fecha y hora actual
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

logging.info(f"Ejecutando vcf_annotator_with_gbff.py versión 6.0 - {current_time}")

# Diccionario para convertir códigos de una letra a tres letras
aa_dict = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe',
    'G': 'Gly', 'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu',
    'M': 'Met', 'N': 'Asn', 'P': 'Pro', 'Q': 'Gln', 'R': 'Arg',
    'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
    '*': 'Ter'
}

def load_gbff(gbff_file):
    logging.info(f"Cargando archivo GBFF: {gbff_file}")
    genes = {}
    contig_mapping = {}
    try:
        for record in SeqIO.parse(gbff_file, "genbank"):
            gbff_contig_name = record.id
            vcf_contig_name = f"AI2644_v1_c{gbff_contig_name.split('_')[1]}"
            contig_mapping[vcf_contig_name] = gbff_contig_name
            
            for feature in record.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
                    product = feature.qualifiers.get("product", ["Unknown"])[0]
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    sequence = str(feature.location.extract(record.seq))
                    if strand == -1:
                        sequence = str(Seq(sequence).reverse_complement())
                    genes[(gbff_contig_name, start, end)] = {
                        "name": gene_name,
                        "product": product,
                        "strand": strand,
                        "sequence": sequence,
                        "start": start,
                        "end": end
                    }
        logging.info(f"Cargados {len(genes)} genes del archivo GBFF")
        logging.info(f"Mapeo de contigs: {contig_mapping}")
    except Exception as e:
        logging.error(f"Error al cargar el archivo GBFF: {str(e)}")
        raise
    return genes, contig_mapping

def get_variant_effect(gene_info, pos, ref, alt):
    logging.debug("Using get_variant_effect function version 6.1")
    
    sequence = gene_info["sequence"]
    strand = gene_info["strand"]
    start = gene_info["start"]
    end = gene_info["end"]
    
    logging.debug(f"Input: pos={pos}, ref={ref}, alt={alt}, strand={strand}, start={start}, end={end}")
    
    if strand == 1:
        rel_pos = pos - start
    else:
        rel_pos = end - pos
    
    codon_start = (rel_pos // 3) * 3
    codon_pos = rel_pos % 3
    
    if strand == -1:
        codon_start = len(sequence) - codon_start - 3
        codon_pos = 2 - codon_pos
    
    ref_codon = sequence[codon_start:codon_start+3]
    alt_codon = list(ref_codon)
    alt_codon[codon_pos] = alt
    alt_codon = "".join(alt_codon)
    
    if strand == -1:
        ref_codon = str(Seq(ref_codon).reverse_complement())
        alt_codon = str(Seq(alt_codon).reverse_complement())
    
    reading_frame = (rel_pos % 3) + 1
    strand_symbol = '+' if strand == 1 else '-'
    
    ref_aa = str(Seq(ref_codon).translate())
    alt_aa = str(Seq(alt_codon).translate())
    
    aa_pos = (rel_pos // 3) + 1
    
    logging.debug(f"Codons: ref_codon={ref_codon}, alt_codon={alt_codon}")
    logging.debug(f"Amino acids: ref_aa={ref_aa}, alt_aa={alt_aa}, aa_pos={aa_pos}")
    
    if aa_pos == 1 and ref_aa == 'M' and alt_aa != 'M':
        aa_change = f"M1?"
        aa_change_3 = f"Met1?"
        variant_type = "start_lost"
    elif ref_aa == alt_aa:
        aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
        aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}{aa_dict[alt_aa]}"
        variant_type = "synonymous_variant"
    elif alt_aa == '*':
        aa_change = f"{ref_aa}{aa_pos}*"
        aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}{aa_dict['*']}"
        variant_type = "stop_gained"
    elif ref_aa == '*':
        aa_change = f"*{aa_pos}{alt_aa}"
        aa_change_3 = f"{aa_dict['*']}{aa_pos}{aa_dict[alt_aa]}"
        variant_type = "stop_lost"
    else:
        aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
        aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}{aa_dict[alt_aa]}"
        variant_type = "missense_variant"
    
    if strand == -1:
        aa_pos_adj = aa_pos  # No ajustamos la posición para la cadena negativa
        aa_change = aa_change.replace(str(aa_pos), str(aa_pos_adj))
        aa_change_3 = aa_change_3.replace(str(aa_pos), str(aa_pos_adj))
    
    logging.debug(f"Final output: variant_type={variant_type}, aa_change={aa_change}, aa_change_3={aa_change_3}")
    
    return variant_type, ref_codon, alt_codon, aa_change, aa_change_3, reading_frame, strand_symbol

def get_indel_effect(gene_info, pos, ref, alt):
    logging.debug(f"Entering get_indel_effect: pos={pos}, ref={ref}, alt={alt}")
    sequence = gene_info["sequence"]
    strand = gene_info["strand"]
    start = gene_info["start"]
    end = gene_info["end"]
    
    logging.debug(f"Gene info: strand={strand}, start={start}, end={end}")
    logging.debug(f"Gene sequence: {sequence[:10]}...{sequence[-10:]}")
    
    if strand == 1:
        rel_pos = pos - start
    else:
        rel_pos = end - pos - 1  # Ajustado para la cadena negativa
    
    aa_pos = rel_pos // 3 + 1
    codon_start = (rel_pos // 3) * 3
    
    logging.debug(f"Calculated positions: rel_pos={rel_pos}, aa_pos={aa_pos}, codon_start={codon_start}")
    
    if strand == 1:
        ref_codon = sequence[codon_start:codon_start+3]
    else:
        codon_start = len(sequence) - codon_start - 3
        ref_codon = str(Seq(sequence[codon_start:codon_start+3]).reverse_complement())
    
    ref_aa = str(Seq(ref_codon).translate())
    reading_frame = (rel_pos % 3) + 1
    strand_symbol = '+' if strand == 1 else '-'
    
    logging.debug(f"Initial ref_codon={ref_codon}, ref_aa={ref_aa}, reading_frame={reading_frame}")
    
    if len(ref) < len(alt):  # Inserción
        inserted_bases = alt[len(ref):]
        if strand == -1:
            inserted_bases = str(Seq(inserted_bases).reverse_complement())
        
        logging.debug(f"Insertion: inserted_bases={inserted_bases}")
        
        if len(inserted_bases) % 3 == 0:
            inserted_aa = str(Seq(inserted_bases).translate())
            aa_change = f"{ref_aa}{aa_pos}ins{inserted_aa}"
            aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}ins{''.join([aa_dict[aa] for aa in inserted_aa])}"
            variant_type = "inframe_insertion"
        else:
            aa_change = f"{ref_aa}{aa_pos}fs"
            aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}fs"
            variant_type = "frameshift_variant"
        alt_codon = "N/A"
    else:  # Deleción
        deleted_bases = ref[len(alt):]
        if strand == -1:
            deleted_bases = str(Seq(deleted_bases).reverse_complement())
        
        logging.debug(f"Deletion: deleted_bases={deleted_bases}")
        
        if len(deleted_bases) % 3 == 0:
            deleted_aa = str(Seq(deleted_bases).translate())
            num_deleted_aa = len(deleted_bases) // 3
            
            if strand == -1:
                aa_pos -= num_deleted_aa - 1
                end_aa_pos = aa_pos + num_deleted_aa - 1
                ref_aa = deleted_aa[0]  # Cambio aquí: usamos el primer aminoácido deletado
            else:
                end_aa_pos = aa_pos + num_deleted_aa - 1
                ref_aa = deleted_aa[0]
            
            logging.debug(f"In-frame deletion: deleted_aa={deleted_aa}, aa_pos={aa_pos}, end_aa_pos={end_aa_pos}, ref_aa={ref_aa}")
            
            aa_change = f"{ref_aa}{aa_pos}_{end_aa_pos}del{deleted_aa}"
            aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}_{end_aa_pos}del{''.join([aa_dict[aa] for aa in deleted_aa])}"
            variant_type = "inframe_deletion"
        else:
            aa_change = f"{ref_aa}{aa_pos}fs"
            aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}fs"
            variant_type = "frameshift_variant"
        alt_codon = "N/A"
    
    logging.debug(f"Final output: variant_type={variant_type}, aa_change={aa_change}, aa_change_3={aa_change_3}")
    
    return variant_type, ref_codon, alt_codon, aa_change, aa_change_3, reading_frame, strand_symbol

def get_variant_type(gene_info, pos, ref, alt):
    if gene_info:
        if len(ref) != len(alt):
            return get_indel_effect(gene_info, pos, ref, alt)
        else:
            return get_variant_effect(gene_info, pos, ref, alt)
    else:
        return "intergenic_variant", "", "", "", "", 0, ""

def get_variant_info(record):
    dp = record.INFO.get('DP', None)
    ad = record.format('AD')[0] if 'AD' in record.FORMAT else None
    af = record.INFO.get('AF', None)
    qual = record.QUAL

    ref_depth = ad[0] if ad is not None else None
    alt_depth = ad[1] if ad is not None and len(ad) > 1 else None

    # AF is already a float in the VCF, so we don't need to convert it
    if af is not None:
        af = af[0] if isinstance(af, (list, tuple)) else af

    return dp, ref_depth, alt_depth, af, qual

def process_vcf(vcf_file, gbff_file, output_file):
    logging.info(f"Procesando VCF: {vcf_file}")
    logging.info(f"Usando GBFF: {gbff_file}")
    logging.info(f"Archivo de salida: {output_file}")

    try:
        genes, contig_mapping = load_gbff(gbff_file)
        vcf_reader = VCF(vcf_file)

        variants_processed = 0
        with open(output_file, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            writer.writerow(['CHROM', 'POS', 'REF', 'ALT', 'GENE', 'PRODUCT', 'REGION_TYPE', 'VARIANT_TYPE', 
                             'REF_CODON', 'ALT_CODON', 'AA_CHANGE', 'AA_CHANGE_3', 'READING_FRAME', 'STRAND', 
                             'TOTAL_DEPTH', 'REF_DEPTH', 'ALT_DEPTH', 'ALT_FREQ', 'QUALITY'])

            for record in vcf_reader:
                vcf_chrom = record.CHROM
                gbff_chrom = contig_mapping.get(vcf_chrom, vcf_chrom)
                pos = record.POS
                ref = record.REF
                alt = str(record.ALT[0])
                
                logging.debug(f"Procesando variante: vcf_chrom={vcf_chrom}, gbff_chrom={gbff_chrom}, pos={pos}, ref={ref}, alt={alt}")
                
                gene_info = None
                gene_name = "Unknown"
                product = "Unknown"
                for (rec_id, start, end), info in genes.items():
                    if rec_id == gbff_chrom and start <= pos <= end:
                        gene_info = info
                        gene_name = info["name"]
                        product = info["product"]
                        break
                
                region_type = "CDS" if gene_info else "intergenic"
                variant_type, ref_codon, alt_codon, aa_change, aa_change_3, reading_frame, strand = get_variant_type(gene_info, pos, ref, alt)
                dp, ref_depth, alt_depth, af, qual = get_variant_info(record)

                # Format the numeric values to ensure consistency
                dp = f"{dp}" if dp is not None else "N/A"
                ref_depth = f"{ref_depth}" if ref_depth is not None else "N/A"
                alt_depth = f"{alt_depth}" if alt_depth is not None else "N/A"
                af = f"{af}" if af is not None else "N/A"  # Simply use the value as is
                qual = f"{qual:.2f}" if qual is not None else "N/A"

                writer.writerow([vcf_chrom, pos, ref, alt, gene_name, product, region_type, variant_type,
                                 ref_codon, alt_codon, aa_change, aa_change_3, reading_frame, strand,
                                 dp, ref_depth, alt_depth, af, qual])
                
                variants_processed += 1
                if variants_processed % 1000 == 0:
                    logging.info(f"Procesadas {variants_processed} variantes")

        logging.info(f"Procesamiento completado. Total de variantes procesadas: {variants_processed}")
        logging.info(f"Resultado guardado en: {output_file}")
    except Exception as e:
        logging.error(f"Error durante el procesamiento: {str(e)}")
        raise

if __name__ == "__main__":
    try:
        if len(sys.argv) != 4:
            raise ValueError("Número incorrecto de argumentos")
        
        vcf_file = sys.argv[1]
        gbff_file = sys.argv[2]
        output_file = sys.argv[3]
        
        if not os.path.exists(vcf_file):
            raise FileNotFoundError(f"El archivo VCF no existe: {vcf_file}")
        if not os.path.exists(gbff_file):
            raise FileNotFoundError(f"El archivo GBFF no existe: {gbff_file}")
        
        process_vcf(vcf_file, gbff_file, output_file)
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        sys.exit(1)