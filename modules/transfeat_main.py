"""
Created on Wed May 15 14:25:00 2019
@author: Juan C Entizne
@email: e.entizne[at]dundee.ac.uk
Modified by: Mojtaba Bagherian
@email: mojtaba.bagheria[at]uwa.edu.au
"""
import os
import sys
import time
import json
import warnings
import logging
from typing import Dict, List, Tuple, Optional, Set, Union

from lib.findlorf.findlorf_tools import *
from lib.transfeat.identify_coding_features import *
from lib.transfeat.identify_non_coding_features import *

from lib.parsing.gtf_object_tools import create_gtf_object, find_utr_regions
from lib.parsing.fasta_parsing_tools import get_fasta_sequences, write_fasta_file
from lib.report.transfeat_report import generate_transfeat_summary

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

warnings.filterwarnings("ignore")

def calculate_exon_only_distance(start_pos: int, end_pos: int, exon_list: List[Tuple[int, int]]) -> int:
    """
    Calculate the distance between two positions excluding introns.
    
    Args:
        start_pos: Start position
        end_pos: End position
        exon_list: List of exon coordinates
        
    Returns:
        int: Distance in nucleotides
    """
    if start_pos > end_pos:
        start_pos, end_pos = end_pos, start_pos

    exon_distance = 0
    for exon in exon_list:
        exon_start, exon_end = min(exon), max(exon)
        if not (exon_end < start_pos or exon_start > end_pos):
            overlap_start = max(exon_start, start_pos)
            overlap_end = min(exon_end, end_pos)
            exon_distance += (overlap_end - overlap_start + 1)
    
    return exon_distance

def calculate_utr_lengths(gtf_obj) -> Dict[str, int]:
    """
    Calculate 3' UTR lengths for all transcripts.
    
    Args:
        gtf_obj: GTF object containing transcript information
        
    Returns:
        Dict[str, int]: Dictionary mapping transcript IDs to their 3' UTR lengths
    """
    threeprimeUTR_len_dt = {}
    
    for trans_id in gtf_obj.trans_exons_dt.keys():
        trans_exons = gtf_obj.trans_exons_dt[trans_id]
        trans_sense = gtf_obj.trans_sense_dt[trans_id]
        trans_cds = gtf_obj.trans_cds_dt.get(trans_id, [])
        
        trans_5utr, trans_3utr, start_codon, stop_codon = find_utr_regions(
            trans_id, trans_sense, trans_exons, trans_cds
        )
        
        if trans_3utr:
            utr3_len = sum([max(exon)-min(exon)+1 for exon in trans_3utr])
            threeprimeUTR_len_dt[trans_id] = utr3_len
        else:
            threeprimeUTR_len_dt[trans_id] = 0
            
    return threeprimeUTR_len_dt

def calculate_junction_distances(
    gtf_obj,
    trans_id: str,
    trans_cds: List[Tuple[int, int]],
    trans_exons: List[Tuple[int, int]],
    trans_sense: str
) -> Tuple[Union[int, str], Union[int, str]]:
    """
    Calculate DSJ and USJ distances for a transcript.
    
    Args:
        gtf_obj: GTF object containing transcript information
        trans_id: Transcript ID
        trans_cds: List of CDS coordinates
        trans_exons: List of exon coordinates
        trans_sense: Transcript strand ('+' or '-')
        
    Returns:
        Tuple[Union[int, str], Union[int, str]]: DSJ and USJ distances
    """
    dsj_distance = "NA"
    usj_distance = "NA"
    
    if not trans_cds:
        return dsj_distance, usj_distance
        
    flat_cds = [c for cds_pair in trans_cds for c in cds_pair]
    if trans_sense == '+':
        trans_start = min(flat_cds)
        trans_stop = max(flat_cds)
    else:
        trans_start = max(flat_cds)
        trans_stop = min(flat_cds)
    
    # Get UTR regions
    trans_5utr, trans_3utr, start_codon, stop_codon = find_utr_regions(
        trans_id, trans_sense, trans_exons, trans_cds
    )
    
    # Calculate USJ distance
    try:
        trans_introns = gtf_obj.trans_introns_dt.get(trans_id, [])
        if trans_introns:
            upstream_junctions = []
            for intron in trans_introns:
                junction_pos = intron[1] if trans_sense == '+' else intron[0]
                is_upstream = (trans_sense == '+' and junction_pos < trans_stop) or \
                             (trans_sense == '-' and junction_pos > trans_stop)
                
                if is_upstream:
                    exon_distance = calculate_exon_only_distance(junction_pos, trans_stop, trans_exons)
                    upstream_junctions.append((junction_pos, exon_distance))
            
            if upstream_junctions:
                closest_junction = min(upstream_junctions, key=lambda x: x[1])
                usj_distance = closest_junction[1]
    except Exception as e:
        logger.error(f"Error calculating UpstreamEJ for {trans_id}: {str(e)}")
    
    # Calculate DSJ distance
    try:
        if trans_3utr and len(trans_3utr) > 1:
            get_introns = lambda exons: [(ex1[-1]+1, ex2[0]-1) for (ex1, ex2) in zip(exons[:-1], exons[1:])]
            introns_3utr = get_introns(trans_3utr)
            
            if introns_3utr:
                # Find the junction closest to the stop codon
                closest_junction = min(
                    introns_3utr,
                    key=lambda x: calculate_exon_only_distance(x[0], trans_stop, trans_exons)
                )
                dsj_distance = calculate_exon_only_distance(closest_junction[0], trans_stop, trans_exons)
    except Exception as e:
        logger.error(f"Error calculating DownstreamEJ for {trans_id}: {str(e)}")
    
    return dsj_distance, usj_distance

def transfeat_main(
    gtf: str,
    fasta: str,
    outpath: str,
    outname: str,
    pep_len: int = 50,
    ptc_len: int = 70,
    uorf_len: int = 10,
    sj_dist: int = 50,
    utr3_len: int = 350,
    orf_index: Optional[str] = None
) -> str:
    """
    Main function for TransFeat analysis with enhanced features.
    
    Args:
        gtf: Path to GTF file
        fasta: Path to FASTA file
        outpath: Output directory path
        outname: Output file name
        pep_len: Minimum peptide length threshold
        ptc_len: PTC length threshold
        uorf_len: uORF length threshold
        sj_dist: Splice junction distance threshold
        utr3_len: 3' UTR length threshold
        orf_index: Path to ORF index file (optional)
        
    Returns:
        str: Path to the generated TransFeat table
    """
    logger.info("Starting TransFeat analysis")

    # +1 AA to account for stop codons during the AA length check
    pep_len += 1

    # In case the user pass the name with a file extension, remove it
    if outname.endswith(".gtf"):
        outname = outname.replace(".gtf", "")
    outname += "_transfeat"

    # Create output folder
    outfolder = os.path.join(outpath, outname)
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    # Get transcriptome annotation
    gtf_obj = create_gtf_object(gtf)

    # Upload transcripts sequences from fasta file
    trans_seq_dt = get_fasta_sequences(fasta)

    # Generate transcripts sequence information
    fasta_header_dt, cds_seq_dt, pep_seq_dt = translate_transcript_cds(trans_seq_dt, gtf_obj)

    logger.info("Retrieving ORF information")
    if orf_index:
        logger.info("Uploading ORF information from ORF index file")
        with open(orf_index) as orf_index_fh:
            orf_dt = json.load(orf_index_fh)
    else:
        # Generate ORF index file
        _, orf_index = find_transcripts_orf_information(gtf, trans_seq_dt, gtf_obj, outfolder)

        logger.info("Uploading ORF information from ORF index file")
        with open(orf_index) as orf_index_fh:
            orf_dt = json.load(orf_index_fh)

    if not orf_dt:
        sys.exit("No ORF information found.")

    # Get transcript start-codon relative position
    relative_start_dt = get_transcript_start_codon_relative_position(gtf_obj)

    # Select authentic stop-codon (at gene level)
    auth_stop_dt = get_genes_authentic_stop_codon_position(gtf_obj)

    logger.info("Retrieving alternative ORFs information")

    # Identify transcripts with long downstream ORF
    is_longer_dorf_dt, ldorf_data_dt = identify_longer_dorf(gtf_obj, relative_start_dt, orf_dt, trans_seq_dt)

    # Identify transcripts with upstream ORF
    is_uorf_dt, uorf_data_dt, urof_categories = identify_uorf(gtf_obj, relative_start_dt, orf_dt, trans_seq_dt, uorf_len)

    logger.info("Identifying Non-Coding features")

    # Identify transcripts without an annotated CDS
    is_orf_absent_dt = is_orf_absent(gtf_obj)

    # Identify transcripts with "Premature Termination Codons" (PTC)
    is_ptc_dt = is_ptc(gtf_obj, ptc_len)

    # Identify transcripts coding for short peptides
    is_orf_short_dt = is_orf_short(gtf_obj, pep_len)

    is_long_3utr_dt = is_long_3utr(gtf_obj, utr3_len)

    # Get transcripts groups for NMD classification
    ptc_trans = set([t_id for t_id, t_bool in is_ptc_dt.items() if t_bool is True])
    long_3utr_trans = set([t_id for t_id, t_bool in is_long_3utr_dt.items() if t_bool is True])
    ov_uorf_trans = urof_categories["overlapping"]
    uorf_trans = urof_categories["not_overlapping"]

    is_nmd_dt, is_dssj_dt = is_nmd(gtf_obj, auth_stop_dt,
                                   sj_dist_th=sj_dist, ptc_trans=ptc_trans,
                                   long_3utr_trans=long_3utr_trans, ov_uorf_trans=ov_uorf_trans, uorf_trans=uorf_trans)

    nmd_features_dt = generate_nmd_features_lines(gtf_obj, is_nmd_dt, is_ptc_dt, is_dssj_dt, is_long_3utr_dt, urof_categories)

    # Identify AS in UTR and NAGNAG features
    as_in_utr_dt, as_utr_location_dt, nagnag_dt = identify_similar_coding_features(gtf_obj)

    # Dictionary of features to annotate
    feature_dicts = {
        "Auto": gtf_obj.trans_gene_dt,
        "No_ORF": is_orf_absent_dt,
        "Short_ORF": is_orf_short_dt,
        "Long_3UTR": is_long_3utr_dt,
        "PTC": is_ptc_dt,
        "NMD": is_nmd_dt,
        "ds_SJ": is_dssj_dt,
        "NMD_features": nmd_features_dt,
        "uORF": urof_categories,
        "ldORF": is_longer_dorf_dt,
        "AS_in_UTR": as_in_utr_dt,
        "AS_Location": as_utr_location_dt,
        "NAGNAG": nagnag_dt
    }

    # Generate the features to annotate into the output table
    coding_potentiality_dt, coding_features_dt, alternative_ORF_dt = generate_feature_tag(gtf_obj, feature_dicts)

    # These dictionaries are required to write the TransFeat table
    table_info_dicts = {
        "Coding_potentiality": coding_potentiality_dt,
        "Coding_features": coding_features_dt,
        "NMD_features": nmd_features_dt,
        "Alternative_ORF": alternative_ORF_dt
    }

    # Get alternative ORFs IDs to validate classification on table
    ldorf_trans = set([t_id for t_id in is_longer_dorf_dt if is_longer_dorf_dt[t_id] is True])
    uorf_trans = set([t_id for t_id in is_uorf_dt if is_uorf_dt[t_id] is True])

    # Calculate 3' UTR lengths and junction distances
    threeprimeUTR_len_dt = calculate_utr_lengths(gtf_obj)
    dsj_distance_dt = {}
    usj_distance_dt = {}

    for trans_id, trans_cds in gtf_obj.trans_cds_dt.items():
        trans_exons = gtf_obj.trans_exons_dt[trans_id]
        trans_sense = gtf_obj.trans_sense_dt[trans_id]
        
        dsj_distance, usj_distance = calculate_junction_distances(
            gtf_obj, trans_id, trans_cds, trans_exons, trans_sense
        )
        
        dsj_distance_dt[trans_id] = dsj_distance
        usj_distance_dt[trans_id] = usj_distance

    # Write detailed splice junction data to separate file with the new column names
    splice_junctions_file = os.path.join(outfolder, f"{outname}_splice_junctions.csv")
    with open(splice_junctions_file, "w+") as fh:
        fh.write("T_ID,Strand,UpstreamEJ,DownstreamEJ,3UTRlength,PTC_dEJ\n")
        for trans_id in sorted(gtf_obj.trans_exons_dt.keys()):
            strand = gtf_obj.trans_sense_dt[trans_id]
            upstream_ej = usj_distance_dt.get(trans_id, "NA")
            downstream_ej = dsj_distance_dt.get(trans_id, "NA")
            utr3_length = threeprimeUTR_len_dt.get(trans_id, 0)
            
            # Determine PTC status based on DownstreamEJ distance
            ptc_dEJ = "No"
            if downstream_ej != "NA":
                try:
                    if int(downstream_ej) >= 50:
                        ptc_dEJ = "Yes"
                except (ValueError, TypeError):
                    pass
            
            fh.write(f"{trans_id},{strand},{upstream_ej},{downstream_ej},{utr3_length},{ptc_dEJ}\n")

    # Write TransFeat table output
    transfeat_table = write_transfeat_table(
        gtf_obj, table_info_dicts, pep_seq_dt, outfolder, outname,
        ldorf_ids=ldorf_trans, uorf_ids=uorf_trans, pep_len=pep_len
    )

    # Write GENERAL/TOTAL output fasta files
    fasta_outfile = os.path.join(outfolder, outname)
    write_fasta_file(cds_seq_dt, f'{fasta_outfile}_nuc.fasta', fasta_header_dt)
    write_fasta_file(pep_seq_dt, f'{fasta_outfile}_pep.fasta', fasta_header_dt)

    # These dictionaries are required to write the fasta files
    sequences_dicts = {
        "Headers": fasta_header_dt,
        "Exonic_seq": trans_seq_dt,
        "CDS_seq": cds_seq_dt,
        "Peptide_seq": pep_seq_dt,
        "ldORF_data": ldorf_data_dt,
        "uORF_data": uorf_data_dt
    }

    # Write ADDITIONAL output fasta files
    write_subcategories_fasta(gtf_obj, transfeat_table, sequences_dicts, outfolder, outname)

    # Generate tables summarizing the main transfeat data
    generate_transfeat_summary(gtf, transfeat_table)

    # Generate table with additional NMD information
    write_NMD_table(gtf_obj, feature_dicts, sequences_dicts, outfolder, outname)

    return transfeat_table
