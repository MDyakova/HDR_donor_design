"""
Functions for working with input data
"""
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
import urllib.request

promoters = pd.read_table('src/static/data/human_epdnew_V33OJ.bed', 
                          header=None, 
                          names=('Chr', 'start', 'end', 'promoter_id', 'v', 'strand'))

def ensemble_info(gene_name):
    """
    Get information from ensemble database.
    Release 109 uses human reference genome GRCh38
    """

    data = EnsemblRelease(109)

    # get gene information
    gene = data.genes_by_name(gene_name)[0]
    strand = gene.strand
    chromosome_name = gene.contig

    gene_dict = gene.to_dict()

    start_gene = gene.start - 1500
    end_gene = gene.end + 1500

    # get sequence
    url = (
        f'https://rest.ensembl.org/sequence/region/human/'
        f'{chromosome_name}:{start_gene}-{end_gene}'
        '?content-type=application/json;version=109'
    )
    data = requests.get(url)
    gene_seq = data.json()["seq"]

    coord_list = list([i for i in range(start_gene, end_gene + 1)])


    # prepare sequence for reverse strand gene location
    if strand == "-":
        compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
        gene_seq = [compl_dict[l] for l in gene_seq][::-1]
        coord_list = coord_list[::-1]

    ensemble_gene_seq = "".join(gene_seq)

    del gene_dict["gene_name"]
    del gene_dict["gene_id"]
    del gene_dict["genome"]

    return ensemble_gene_seq, gene_dict, strand


def ncbi_information(ncbi_id):
    """
    Get information from NCBI database
    """

    transcripts_info = []
    Entrez.email = "Your.Name.Here@example.org"

    # get information about transcript
    handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="gb", retmode="text")
    seq_record = [seq_record for seq_record in SeqIO.parse(handle, "gb")][0]
    refseq_sequence = str(seq_record.seq)
    features_list = []
    for feature in seq_record.features:
        features_list.append(
            [
                feature.type,
                int(feature.location.start + 1),
                int(feature.location.end),
                feature.qualifiers,
            ]
        )
        if feature.type == "CDS":
            protein_id = feature.qualifiers["protein_id"][0]

    annotation = seq_record.annotations["keywords"]
    ensemble_is_rs = ""
    if "MANE Select" in annotation:
        ensemble_is_rs = (
            seq_record.annotations["structured_comment"]["RefSeq-Attributes"][
                "MANE Ensembl match"
            ]
            .split("/")[0]
            .split(".")[0]
        )

    transcripts_info.append(
        [ncbi_id, refseq_sequence, features_list, ensemble_is_rs, protein_id]
    )
    handle.close()

    # get information about coding sequence
    for rna_info in transcripts_info:
        features = rna_info[2]
        ensemble_is_rs = rna_info[3]
        protein_id = rna_info[4]
        cds = list(filter(lambda p: p[0] == "CDS", features))[0]
        cds_start_pos = cds[1]
        cds_end_pos = cds[2]

    cds_seq = refseq_sequence[cds_start_pos - 1 : cds_end_pos]
    cds_seq_long = refseq_sequence[cds_start_pos - 22 : cds_end_pos + 21]

    return transcripts_info, cds_seq, cds_seq_long


def guide_info(guide_seq, cds_seq, strand, ensemble_gene_seq, cds_seq_long):
    """
    Search guide position and start codon in gene.
    """

    compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    if ';' not in guide_seq:
        if strand == '-':
            guide = ''.join([compl_dict[i] for i in guide_seq][::-1])
        else:
            guide = guide_seq
    else:
        left_guide, right_guide = guide_seq.split(';')
        left_guide = left_guide.strip()
        right_guide = right_guide.strip()
        
        if left_guide not in ensemble_gene_seq:
            left_guide = ''.join([compl_dict[i] for i in left_guide][::-1])
        if right_guide not in ensemble_gene_seq:
            right_guide = ''.join([compl_dict[i] for i in right_guide][::-1])
        
        left_guide_start = len(ensemble_gene_seq.split(left_guide)[0])
        right_guide_start = len(ensemble_gene_seq.split(right_guide)[0])
        
        if left_guide_start>right_guide_start:
            guide = (right_guide 
                    + ensemble_gene_seq.split(right_guide)[1].split(left_guide)[0] 
                    + left_guide)
        else:
            guide = (left_guide 
                    + ensemble_gene_seq.split(left_guide)[1].split(right_guide)[0] 
                    + right_guide)  

    # search cut site
    guide_cut_size = len(guide) // 2

    # search position in reading frame
    if guide[:guide_cut_size] in cds_seq:
        guide_in_cds_pos = len(cds_seq.split(guide[:guide_cut_size])[0])
        guide_pos = 0
    elif guide[-guide_cut_size:] in cds_seq:
        guide_in_cds_pos = len(cds_seq.split(guide[-guide_cut_size:])[0])
        guide_pos = len(guide.split(guide[-guide_cut_size:])[0])
    else:
        print("guide not found")
        # CDS_seq = refseq_sequence[CDS_start_pos-1-21:CDS_end_pos+21]
        if guide[:guide_cut_size] in cds_seq_long:
            # guide_in_cds_seq = guide[:guide_cut_size]
            guide_in_cds_pos = len(cds_seq_long.split(guide[:guide_cut_size])[0])
            guide_pos = 0
        elif guide[-guide_cut_size:] in cds_seq_long:
            # guide_in_cds_seq = guide[-guide_cut_size:]
            guide_in_cds_pos = len(cds_seq_long.split(guide[-guide_cut_size:])[0])
            guide_pos = len(guide.split(guide[-guide_cut_size:])[0])

    codon_positions = {}
    for i, k in zip(range(len(guide)), range(len(guide))):
        codon_positions[k] = (guide_in_cds_pos - (guide_pos - i)) % 3

    # cut site position in frame
    cut_site_codon_pos = codon_positions[guide_cut_size]
    position_insert_start = guide_cut_size - cut_site_codon_pos

    return position_insert_start, guide_cut_size, guide, guide_in_cds_pos


def make_seqience(
    flank_size, 
    lha_size, 
    rha_size, 
    donor_elements, 
    element_sequnces, 
    colors, 
    guide_in_cds_pos, 
    promoter_list,
    left_seq,
    cds_seq
):
    """
    Make insert sequence for all selected elements.
    Make colored sequence.
    Search new start codon if necessary.
    """
    stop_codons = ["TAA", "TAG", "TGA"]
    compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    insert_start = flank_size + lha_size + 1
    next_element = insert_start

    elements_list = []
    elements_list.append(
        ["LHA", flank_size + 1, flank_size + lha_size, "+", "gene sequence"]
    )

    atg_20_seq = cds_seq[cds_seq.index('ATG'):][:20]
    if atg_20_seq in left_seq:
        atg_start = len(left_seq.split(atg_20_seq)[0])
        if atg_start<insert_start:
            elements_list.append(['ATG_gene', atg_start, atg_start + 3, '+', 'Start codon'])

     
    for p_sequence in promoter_list:
        if p_sequence in left_seq:
            promoter_start = len(left_seq.split(p_sequence)[0])
            promoter_end = promoter_start + len(p_sequence)
            elements_list.append(['Promoter_gene', promoter_start, promoter_end, '+', 'Promoter'])

        else:
            for delta in range(1, len(p_sequence)//2):
                if p_sequence[delta:] in left_seq:
                    promoter_start = len(left_seq.split(p_sequence[delta:])[0])
                    promoter_end = promoter_start + len(p_sequence[delta:])
                    elements_list.append(['Promoter_gene', promoter_start, promoter_end, '+', 'Promoter'])

    insert_sequence = ""
    insert_sequence_color = ""
    coding_sequence = ""
    new_start_codon = False
    for element in donor_elements:
        element = "_".join(element.split("_")[1:])
        if "_reverse" in element:
            direction = "-"
            element = element.split("_reverse")[0]
            seq_i = (
                element_sequnces[
                    element_sequnces["Names"] == element
                ]["Sequence"]
                .max()
                .upper()
            )
            group = element_sequnces[
                element_sequnces["Names"] == element
            ]["Elements"].max()
            seq_i = seq_i.replace("U", "T")
            seq_i = "".join([compl_dict[i] for i in seq_i][::-1])
        else:
            direction = "+"
            seq_i = (
                element_sequnces[
                    element_sequnces["Names"] == element
                ]["Sequence"]
                .max()
                .upper()
            )
            group = element_sequnces[
                element_sequnces["Names"] == element
            ]["Elements"].max()
            seq_i = seq_i.replace("U", "T")

        in_frame = element_sequnces[
            element_sequnces["Names"] == element
        ]["in frame"].max()
        if (seq_i[-3:] in stop_codons) & (in_frame == 1):
            seq_i = seq_i[:-3]
            print(element, "STOP")

        insert_sequence += seq_i

        if elements_list[-1][4] == "Promoter":
            coding_sequence = ""
            new_start_codon = True
        coding_sequence += seq_i

        insert_sequence_color += colors[group][2]
        insert_sequence_color += seq_i
        insert_sequence_color += "</span>"

        if (len(coding_sequence) % 3 != 0) & (in_frame == 1):
            if new_start_codon:
                atg_ind = np.argmax(
                    [
                        coding_sequence[i : i + 3] == "ATG"
                        for i in range(len(coding_sequence))
                    ]
                )
                new_start_codon = False
                if len(coding_sequence[atg_ind:]) % 3 != 0:
                    n_tail = (3 - len(coding_sequence[atg_ind:]) % 3)
                    insert_sequence += "N" * n_tail
                    coding_sequence += "N" * n_tail
                    insert_sequence_color += colors[group][2]
                    insert_sequence_color += "N" * n_tail
                    insert_sequence_color += "</span>"
            else:
                n_tail = (3 - len(insert_sequence) % 3)
                insert_sequence += "N" * n_tail
                insert_sequence_color += colors[group][2]
                insert_sequence_color += "N" * n_tail
                insert_sequence_color += "</span>"

        elements_list.append(
            [
                element,
                next_element,
                insert_start + len(insert_sequence) - 1,
                direction,
                group,
            ]
        )
        next_element = insert_start + len(insert_sequence)

    elements_list.append(
        [
            "RHA",
            elements_list[-1][2] + 1,
            elements_list[-1][2] + rha_size,
            "+",
            "gene sequence",
        ]
    )

    return elements_list, insert_sequence, insert_sequence_color


def make_sequence_image(gene_name, elements_list, colors, full_sequence):
    """
    Make image of all elements.
    """
    features = []
    plt.figure(figsize=(10, 20))
    for element in elements_list:
        if element[3] == "+":
            strand_plot = +1
        else:
            strand_plot = -1
        feature = GraphicFeature(
            start=element[1],
            end=element[2],
            strand=strand_plot,
            color=colors[element[4]][0],
            label=element[0],
        )
        features.append(feature)
    #     break
    record = GraphicRecord(sequence_length=len(full_sequence), features=features)
    cropped_record = record.crop((1, len(full_sequence)))
    cropped_record.sequence = full_sequence
    cropped_record.plot(
        figure_width=10, strand_in_label_threshold=7, plot_sequence=False
    )
    plt.savefig(
        "src/static/outputs/" + gene_name + "/map_for_" + gene_name + ".png",
        bbox_inches="tight",
    )


def save_files(
    gene_name,
    elements_list,
    full_sequence,
    colors,
    files_name
):
    """
    Save all files for SnapGene
    """
    date_today = str(date.today())
    features_list_donor = []

    for feature in elements_list:
        start = feature[1]
        end = feature[2]
        name = feature[0]
        color = colors[feature[4]][0]
        direction = feature[3]
        features_list_donor.append(
            [
                "0",
                int(start),
                int(end),
                name,
                1000,
                direction,
                0,
                1,
                ",".join([str(int(c * 255)) for c in m_colors.to_rgb(color)]),
                None,
                None,
                None,
            ]
        )

    features_df = pd.DataFrame(features_list_donor)
    features_df[1] = features_df[1] - 1

    features_df.to_csv(
        "src/static/outputs/" + gene_name + "/" + files_name + "_donor.bed",
        sep="\t",
        header=None,
        index=None,
    )

    with open("src/static/outputs/" + gene_name + "/" + files_name + "_donor.bed",
              "r", encoding="utf-8") as file:
        bed_file = file.read()

    bed_file_name = (
        "src/static/outputs/"
        + gene_name
        + "/"
        + files_name
        + ".bed"
    )
    with open(bed_file_name, "w", encoding="utf-8") as file:
        file.write("track itemRgb=On\n")
        file.write(bed_file)

    fasta_file_name = (
        "src/static/outputs/" + gene_name + "/" + files_name + "_donor_sequence.fa"
    )
    with open(fasta_file_name, "w", encoding="utf-8") as file:
        file.write("> " + gene_name + "\n")
        file.write(full_sequence + "\n")

    return fasta_file_name, bed_file_name

def list_files_and_sizes(folder_path):
    files_and_sizes = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)
            file_size = str(int(os.path.getsize(file_path) / 1024.0)) + 'kb'
            files_and_sizes.append((file_path, file_size))

    return files_and_sizes

def find_promoter(gene_name):
    """
    Find promoter sequence in LHA
    """
    promoter_gene = promoters[promoters['promoter_id'].apply(lambda p: gene_name == p.split('_')[0])]

    promoter_list = [] 
    p_sequence = ''
    for promoter_id in pd.unique(promoter_gene['promoter_id']):
        url = 'https://epd.expasy.org/cgi-bin/epd/get_doc?db=hgEpdNew&format=genome&entry=' + promoter_id

        fp = urllib.request.urlopen(url)
        mybytes = fp.read()
        html = mybytes.decode("utf8")
        fp.close()

        p_sequence = html.split('Sequence:</td><td>')[1].split('</td></tr>')[0]
        p_sequence = p_sequence.upper()
        
        promoter_list.append(p_sequence)

    return promoter_list
