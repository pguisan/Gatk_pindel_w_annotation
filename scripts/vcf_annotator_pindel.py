#!/usr/bin/env python3
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from cyvcf2 import VCF
from datetime import datetime
import logging

# Configurar logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Obtener la fecha y hora actual
current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
logging.info(f"Ejecutando vcf_annotator_pindel.py versión 1.2 - {current_time}")

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
                genes[(gbff_contig_name, start, end)] = {
                    "name": gene_name,
                    "product": product,
                    "strand": strand,
                    "sequence": sequence,
                    "start": start,
                    "end": end
                }
    logging.info(f"Cargados {len(genes)} genes del archivo GBFF")
    return genes, contig_mapping

def translate_sequence(sequence, strand):
    if strand == -1:
        sequence = str(Seq(sequence).reverse_complement())
    return str(Seq(sequence).translate())

def get_variant_effect(gene_info, pos, ref, alt):
    sequence = gene_info["sequence"]
    strand = gene_info["strand"]
    start = gene_info["start"]
    end = gene_info["end"]
    
    if strand == 1:
        rel_pos = pos - start
    else:
        rel_pos = end - pos - 1  # Ajustado para la cadena negativa
    
    aa_pos = rel_pos // 3 + 1
    codon_start = (rel_pos // 3) * 3
    
    if strand == 1:
        ref_codon = sequence[codon_start:codon_start+3]
    else:
        codon_start = len(sequence) - codon_start - 3
        ref_codon = str(Seq(sequence[codon_start:codon_start+3]).reverse_complement())
    
    ref_aa = str(Seq(ref_codon).translate())
    reading_frame = (rel_pos % 3) + 1
    strand_symbol = '+' if strand == 1 else '-'
    
    if len(ref) < len(alt):  # Inserción
        inserted_seq = alt[len(ref):]
        if strand == -1:
            inserted_seq = str(Seq(inserted_seq).reverse_complement())
        
        if len(inserted_seq) % 3 == 0:
            inserted_aa = translate_sequence(inserted_seq, 1)
            aa_change = f"{ref_aa}{aa_pos}ins{inserted_aa}"
            aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}ins{''.join([aa_dict[aa] for aa in inserted_aa])}"
            variant_type = "inframe_insertion"
        else:
            aa_change = f"{ref_aa}{aa_pos}fs"
            aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}fs"
            variant_type = "frameshift_variant"
        alt_codon = "N/A"
    elif len(ref) > len(alt):  # Deleción
        deleted_seq = ref[len(alt):]
        if strand == -1:
            deleted_seq = str(Seq(deleted_seq).reverse_complement())
        
        if len(deleted_seq) % 3 == 0:
            deleted_aa = translate_sequence(deleted_seq, 1)
            if len(deleted_aa) == 1:
                aa_change = f"del{deleted_aa}{aa_pos}"
                aa_change_3 = f"del{aa_dict[deleted_aa]}{aa_pos}"
            else:
                end_pos = aa_pos + len(deleted_aa) - 1
                aa_change = f"del{deleted_aa[0]}{aa_pos}_{deleted_aa[-1]}{end_pos}"
                aa_change_3 = f"del{aa_dict[deleted_aa[0]]}{aa_pos}_{aa_dict[deleted_aa[-1]]}{end_pos}"
            variant_type = "inframe_deletion"
        else:
            aa_change = f"{ref_aa}{aa_pos}fs"
            aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}fs"
            variant_type = "frameshift_variant"
        alt_codon = "N/A"
    else:  # SNP (aunque es poco probable en resultados de Pindel)
        alt_codon = list(ref_codon)
        alt_codon[rel_pos % 3] = alt
        alt_codon = "".join(alt_codon)
        
        if strand == -1:
            alt_codon = str(Seq(alt_codon).reverse_complement())
        
        alt_aa = translate_sequence(alt_codon, 1)
        
        aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
        aa_change_3 = f"{aa_dict[ref_aa]}{aa_pos}{aa_dict[alt_aa]}"
        
        if ref_aa == alt_aa:
            variant_type = "synonymous_variant"
        elif alt_aa == '*':
            variant_type = "stop_gained"
        elif ref_aa == '*':
            variant_type = "stop_lost"
        else:
            variant_type = "missense_variant"
    
    logging.debug(f"Variant effect: pos={pos}, ref={ref}, alt={alt}, ref_codon={ref_codon}, aa_pos={aa_pos}, reading_frame={reading_frame}, strand={strand}")
    
    return variant_type, ref_codon, alt_codon, aa_change, aa_change_3, reading_frame, strand_symbol

def find_affected_gene(genes, chrom, pos):
    """
    Encuentra el gen que contiene una posición específica.
    Útil para determinar el gen exacto en el que comienza una variante.
    """
    for (rec_id, start, end), info in genes.items():
        if rec_id == chrom and start <= pos <= end:
            return info
    return None

def find_affected_genes(genes, chrom, start, end):
    """
    Encuentra todos los genes afectados por un rango de posiciones.
    Útil para variantes grandes que pueden afectar a múltiples genes.
    """
    affected_genes = []
    for (rec_id, gene_start, gene_end), info in genes.items():
        if rec_id == chrom and (
            (start <= gene_start <= end) or
            (start <= gene_end <= end) or
            (gene_start <= start and end <= gene_end)
        ):
            affected_genes.append((gene_start, info['name']))
    return sorted(affected_genes)

def process_pindel_output(input_file, output_file, gbff_file, min_indel_size=50):
    genes, contig_mapping = load_gbff(gbff_file)
    
    filtered_count = 0
    total_count = 0
    
    with open(output_file, 'w', newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['CHROM', 'POS', 'REF', 'ALT', 'SVLEN', 'AFFECTED_GENES', 'PRODUCTS', 'REGION_TYPE', 'VARIANT_TYPE', 
                         'REF_CODON', 'ALT_CODON', 'AA_CHANGE', 'AA_CHANGE_3', 'READING_FRAME', 'STRAND', 
                         'PINDEL_SVTYPE', 'PINDEL_END', 'SUPPORTING_READS', 'PINDEL_FILTER'])
        
        vcf_reader = VCF(input_file)
        for record in vcf_reader:
            total_count += 1
            vcf_chrom = record.CHROM
            gbff_chrom = contig_mapping.get(vcf_chrom, vcf_chrom)
            pos = record.POS
            ref = record.REF
            alt = str(record.ALT[0])
            
            svlen = abs(len(alt) - len(ref))
            
            if svlen < min_indel_size:
                filtered_count += 1
                continue
            
            pindel_filter = record.FILTER if record.FILTER else "PASS"
            svtype = record.INFO.get('SVTYPE', 'Unknown')
            end = record.INFO.get('END', pos + svlen)
            supporting_reads = record.INFO.get('SR', 'N/A')
            
            affected_genes = find_affected_genes(genes, gbff_chrom, pos, end)
            
            if affected_genes:
                first_gene = affected_genes[0][1]
                affected_gene = find_affected_gene(genes, gbff_chrom, pos)
                if affected_gene:
                    gene_name = affected_gene['name']
                    region_type = "CDS"
                    variant_type, ref_codon, alt_codon, aa_change, aa_change_3, reading_frame, strand = get_variant_effect(affected_gene, pos, ref, alt)
                else:
                    gene_name = "Intergenic"
                    region_type = "intergenic"
                    variant_type, ref_codon, alt_codon, aa_change, aa_change_3, reading_frame, strand = "intergenic_variant", "", "", "", "", 0, ""
                
                affected_genes_str = "; ".join([gene[1] for gene in affected_genes])
                products = "; ".join([genes.get((gbff_chrom, g[0], g[0]), {}).get('product', 'N/A') for g in affected_genes])
            else:
                affected_genes_str = "Intergenic"
                gene_name = "Intergenic"
                products = "N/A"
                region_type = "intergenic"
                variant_type, ref_codon, alt_codon, aa_change, aa_change_3, reading_frame, strand = "intergenic_variant", "", "", "", "", 0, ""

            writer.writerow([vcf_chrom, pos, ref, alt, svlen,
                             affected_genes_str, products, region_type, variant_type,
                             ref_codon, alt_codon, aa_change, aa_change_3, reading_frame, strand,
                             svtype, end, supporting_reads, pindel_filter])

    logging.info(f"Procesamiento completado. Resultado guardado en: {output_file}")
    logging.info(f"Total de variantes procesadas: {total_count}")
    logging.info(f"Variantes filtradas (tamaño < {min_indel_size}): {filtered_count}")
    logging.info(f"Variantes reportadas: {total_count - filtered_count}")

if __name__ == "__main__":
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print("Uso: python vcf_annotator_pindel.py <input_file> <output_file> <gbff_file> [min_indel_size]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    gbff_file = sys.argv[3]
    min_indel_size = int(sys.argv[4]) if len(sys.argv) == 5 else 50
    
    logging.info(f"Iniciando procesamiento con tamaño mínimo de indel: {min_indel_size}")
    process_pindel_output(input_file, output_file, gbff_file, min_indel_size)