import os
from datetime import date
import pandas as pd
import numpy as np
from pyensembl import EnsemblRelease
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
from Bio import Entrez, SeqIO
from matplotlib import colors as m_colors
import requests

from utilities import (
    ensemble_info,
    ncbi_information,
    guide_info,
    make_seqience,
    make_sequence_image,
    save_files,
    list_files_and_sizes
)

# table with all possible elements and sequences
element_sequences = pd.read_excel(
    "data/all_sequences.xlsx", sheet_name="Sequences"
)
element_sequences = element_sequences.groupby(
    by=["Elements", "Names"], as_index=False
).max()
element_sequnces = element_sequences.sort_values(
    by=["Elements", "Names"], key=lambda col: col.str.lower()
)
possible_elements = [
    {"group": i, "name": j}
    for i, j in zip(
        element_sequences["Elements"], element_sequences["Names"]
    )
]

colors = {'2A motif':['#0000EE', '(0, 0, 238)', "<span class='blue-text'>"],
          'protein': ['#00C957', '(0, 201, 87)', "<span class='green-text'>"],
          'cloning' : ['#CDB38B', '(205, 179, 139)', "<span class='peach-text'>"],
          'Stop codon': ['#FF3030', '(255, 48, 48)', "<span class='red-text'>"],
          'Terminator': ['#FF6103', '(255, 97, 3)', "<span class='orange-text'>"],
          'custom':['#00C957', '(0, 201, 87)', "<span class='green-text'>"],
          'Promoter':['#e0441d', '(224, 68, 29)', "<span class='redlight-text'>"],
          'signal peptide':['#8b1de0', '(139, 29, 224)', "<span class='purple-text'>"],
          'CAP binding site':['#8b1de0', '(139, 29, 224)', "<span class='purple-text'>"],
          'Kozak sequence':['#8b1de0', '(139, 29, 224)', "<span class='purple-text'>"],
          'transport':['#8b1de0', '(139, 29, 224)', "<span class='purple-text'>"],
          'gene sequence':['#b5b5b1', '(181, 181, 177)', "<span class='grey-text'>"],
          '5UTR':['#3737c4', '(55, 55, 196)', "<span class='blue-light-text'>"],
          'new':['#050505', '(0, 0, 0)', "<span class='black-text'>"]}

def test_ensemble_info():
    """Check information from ensemble database"""

    gene_name = 'PDCD1'
    first_20_nucleotides = 'GGCTGGGGGCAGAGGGAGGT'

    ensemble_gene_seq, gene_dict, strand = ensemble_info(gene_name)

    assert (
        ensemble_gene_seq[:20] == first_20_nucleotides
    ), "Ensemble data isn't correct"

def test_ncbi_information():
    """Check information from NCBI database"""

    ncbi_id = 'NM_005018.3'
    protein_id = 'NP_005009.2'

    transcripts_info, cds_seq = ncbi_information(ncbi_id)

    assert (
        transcripts_info[0][-1] == protein_id
    ), "NCBI data isn't correct"

def test_guide_info():
    """Check guide information"""
    guide_seq = 'AGTTGTAGCACCGCCCAGACGACTGGCCAGGGCGCCTGTGGGATCTGCATGCCTG'
    gene_name = 'PDCD1'
    ncbi_id = 'NM_005018.3'
    start_position = 26
    cut_size = 27 
    guide_seq_correct = 'CAGGCATGCAGATCCCACAGGCGCCCTGGCCAGTCGTCTGGGCGGTGCTACAACT'

    _, _, strand = ensemble_info(gene_name)
    _, cds_seq = ncbi_information(ncbi_id)

    position_insert_start, guide_cut_size, guide = guide_info(guide_seq, cds_seq, strand)

    assert (
        position_insert_start == start_position
    ), "Incorrect start possition"
    assert (
        guide_cut_size == cut_size
    ), "Incorrect cut size"
    assert (
        guide == guide_seq_correct
    ), "Incorrect guide sequence"

def test_make_seqience():
    """Check maked sequence"""

    flank_size = 100
    lha_size = 400
    rha_size = 400
    donor_elements = ['2A motif_T2A']

    check_element = ['T2A', 501, 554, '+', '2A motif']
    check_seq = 'GAGGGCAGAGGCAGTCTGCTGACATGCGGTGACGTGGAAGAGAATCCCGGCCCT'
    check_seq_color = "<span class='blue-text'>GAGGGCAGAGGCAGTCTGCTGACATGCGGTGACGTGGAAGAGAATCCCGGCCCT</span>"

    elements_list, insert_sequence, insert_sequence_color = make_seqience(flank_size, 
                                                                        lha_size, 
                                                                        rha_size, 
                                                                        donor_elements, 
                                                                        element_sequnces, 
                                                                        colors)

    assert (
        elements_list[1] == check_element
    ), "Incorrect element list"

    assert (
        insert_sequence == check_seq
    ), "Incorrect insert sequence"

    assert (
        insert_sequence_color == check_seq_color
    ), "Incorrect insert color sequence"