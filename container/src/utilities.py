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

amino_acid_dict = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
    "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Stop": "*"
}

codon_using = pd.read_table('src/static/data/nuclear_codon_statistics.tsv')
# codon_using = codon_using[codon_using['Fraction']>=0.08]
# codon_using = codon_using[codon_using['Preferred Codon']=="True"]
codon_using['aa'] = codon_using['Amino acid'].apply(lambda p: amino_acid_dict[p])
codon_using.sort_values(by=['aa', 'Fraction'], ascending=False, inplace=True)

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
    exon_number = 1
    for feature in seq_record.features:
        feature_seq = refseq_sequence[feature.location.start:feature.location.end]
        if feature.type == 'exon':
            features_list.append([f'exon_{exon_number}', 
                                  int(feature.location.start+1), 
                                  int(feature.location.end), 
                                  feature.qualifiers, feature_seq])
            exon_number += 1
        else:
            features_list.append([feature.type, 
                                  int(feature.location.start+1), 
                                  int(feature.location.end), 
                                  feature.qualifiers, feature_seq])
        
        if feature.type == 'CDS':
            protein_id = feature.qualifiers['protein_id'][0]

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

    return transcripts_info, cds_seq, cds_seq_long, features_list, refseq_sequence


def guide_info(guide_seq, cds_seq, strand, ensemble_gene_seq, cds_seq_long, refseq_sequence):
    """
    Search guide position and start codon in gene.
    """

    compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}

    if ';' not in guide_seq:
        right_guide = guide_seq[-20:]
        left_guide = ''.join([compl_dict[i] for i in guide_seq[:20]][::-1])
        guide_seq = ';'.join([left_guide, right_guide])

        right_guide_save = right_guide
        left_guide_save = left_guide

    if ';' not in guide_seq:
        if strand == '-':
            guide = ''.join([compl_dict[i] for i in guide_seq][::-1])
        else:
            guide = guide_seq
    else:
        left_guide, right_guide = guide_seq.split(';')
        left_guide = left_guide.strip()
        right_guide = right_guide.strip()

        right_guide_save = right_guide
        left_guide_save = left_guide
        
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
        guide_in_cds_seq = guide[:guide_cut_size]
        guide_in_cds_pos = len(cds_seq.split(guide[:guide_cut_size])[0])
        guide_pos = 0
    elif guide[-guide_cut_size:] in cds_seq:
        guide_in_cds_seq = guide[-guide_cut_size:]
        guide_in_cds_pos = len(cds_seq.split(guide[-guide_cut_size:])[0])
        guide_pos = len(guide.split(guide[-guide_cut_size:])[0])
    else:
        print("guide not found")
        # CDS_seq = refseq_sequence[CDS_start_pos-1-21:CDS_end_pos+21]
        if guide[:guide_cut_size] in cds_seq_long:
            guide_in_cds_seq = guide[:guide_cut_size]
            guide_in_cds_pos = len(cds_seq_long.split(guide[:guide_cut_size])[0])
            guide_pos = 0
        elif guide[-guide_cut_size:] in cds_seq_long:
            guide_in_cds_seq = guide[-guide_cut_size:]
            guide_in_cds_pos = len(cds_seq_long.split(guide[-guide_cut_size:])[0])
            guide_pos = len(guide.split(guide[-guide_cut_size:])[0])

    codon_positions = {}
    for i, k in zip(range(len(guide)), range(len(guide))):
        codon_positions[k] = (guide_in_cds_pos - (guide_pos - i)) % 3

    # cut site position in frame
    cut_site_codon_pos = codon_positions[guide_cut_size]
    position_insert_start = guide_cut_size - cut_site_codon_pos
    guide_in_transcript_pos = len(refseq_sequence.split(guide_in_cds_seq)[0])

    return (position_insert_start, 
            guide_cut_size, guide, 
            guide_in_cds_pos, 
            guide_in_cds_seq, 
            guide_in_transcript_pos, 
            guide_pos, left_guide_save, right_guide_save)


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
    next_element_no_n = insert_start

    elements_list = []
    elements_list_no_n = []
    elements_list.append(
        ["LHA", flank_size + 1, flank_size + lha_size, "+", "gene sequence"]
    )
    elements_list_no_n.append(
        ["LHA", flank_size + 1, flank_size + lha_size, "+", "gene sequence"]
    )

    atg_20_seq = cds_seq[cds_seq.index('ATG'):][:20]
    if atg_20_seq in left_seq:
        atg_start = len(left_seq.split(atg_20_seq)[0])
        if atg_start<insert_start:
            elements_list.append(['ATG_gene', atg_start+1, atg_start + 3, '+', 'Start codon'])
            elements_list_no_n.append(['ATG_gene', atg_start+1, atg_start + 3, '+', 'Start codon'])

     
    for p_sequence in promoter_list:
        if p_sequence in left_seq:
            promoter_start = len(left_seq.split(p_sequence)[0])
            promoter_end = promoter_start + len(p_sequence)
            elements_list.append(['Promoter_gene', promoter_start+1, promoter_end, '+', 'Promoter'])
            elements_list_no_n.append(['Promoter_gene', promoter_start+1, promoter_end, '+', 'Promoter'])

        else:
            for delta in range(1, len(p_sequence)//2):
                if p_sequence[delta:] in left_seq:
                    promoter_start = len(left_seq.split(p_sequence[delta:])[0])
                    promoter_end = promoter_start + len(p_sequence[delta:])
                    elements_list.append(['Promoter_gene', promoter_start+1, promoter_end, '+', 'Promoter'])
                    elements_list_no_n.append(['Promoter_gene', promoter_start+1, promoter_end, '+', 'Promoter'])

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
            ])
        elements_list_no_n.append(
            [
                element,
                next_element_no_n,
                insert_start + len(insert_sequence.replace('N', '')) - 1,
                direction,
                group,
            ])
        next_element = insert_start + len(insert_sequence)
        next_element_no_n = insert_start + len(insert_sequence.replace('N', ''))

    elements_list.append(
        [
            "RHA",
            np.max([i[2] for i in elements_list]) + 1,
            np.max([i[2] for i in elements_list]) + rha_size,
            "+",
            "gene sequence",
        ]
    )
    elements_list_no_n.append(
        [
            "RHA",
            np.max([i[2] for i in elements_list_no_n]) + 1,
            np.max([i[2] for i in elements_list_no_n]) + rha_size,
            "+",
            "gene sequence",
        ]
    )

    return elements_list, insert_sequence, insert_sequence_color, elements_list_no_n


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

def gene_features(features_list, guide, guide_pos, guide_in_transcript_pos, 
                  elements_list, position_insert_start, delta_nucleotides):
    exons = list(filter(lambda p: 'exon' in p[0], features_list))
    not_exons = list(filter(lambda p: 'exon' not in p[0], features_list[2:]))

    new_exons_features = []
    for exon in exons:
        exon_start = exon[1]
        exon_end = exon[2]
        exon_seq = exon[4]
        for feature in not_exons:
            feature_start = feature[1]
            feature_end = feature[2]
            feature_seq = feature[4]
            if feature_end<exon_start:
                pass
            elif feature_start>exon_end:
                pass
            else:
                if feature[0]=='misc_feature':
                    new_feature_start = np.maximum(feature_start, exon_start)
                    new_feature_end = np.minimum(feature_end, exon_end)
                    new_exons_features.append([exon[0], 
                                               feature[3]['note'][0].split(';')[-1].split('/')[0].strip(), 
                                               new_feature_start, new_feature_end, feature[3]])
                else:
                    new_feature_start = np.maximum(feature_start, exon_start)
                    new_feature_end = np.minimum(feature_end, exon_end)
                    new_exons_features.append([exon[0], feature[0], new_feature_start, new_feature_end, feature[3]])



    exon = list(filter(lambda p: (p[1]<=guide_in_transcript_pos) & (p[2]>=guide_in_transcript_pos), exons))[0]
    exon_name = exon[0]
    new_exons_features = list(filter(lambda p: p[0]==exon_name, new_exons_features))

    new_exons_features_before = []
    new_exons_features_after = []
    for feature in new_exons_features:
        feature_start = feature[2]
        feature_end = feature[3]
        if feature_end<(guide_in_transcript_pos + position_insert_start - guide_pos):
            new_exons_features_before.append(feature)
            pass
        elif feature_start>(guide_in_transcript_pos + (len(guide) - position_insert_start - 1) - guide_pos):
            new_exons_features_after.append(feature)
        else:
            new_feature_start = np.maximum(feature_start, 
                                           (guide_in_transcript_pos + position_insert_start - guide_pos) + delta_nucleotides + 1)
            new_feature_end = np.minimum(feature_end, (guide_in_transcript_pos + position_insert_start - guide_pos))
            new_exons_features_before.append([feature[0], feature[1], feature_start, new_feature_end, feature[4]])
            new_exons_features_after.append([feature[0], feature[1], new_feature_start, feature_end, feature[4]])

    new_feature_start = np.maximum(exon[1], (guide_in_transcript_pos + position_insert_start - guide_pos) + delta_nucleotides + 1)
    new_feature_end = np.minimum(exon[2], (guide_in_transcript_pos + position_insert_start - guide_pos))
    new_exons_features_before.append([exon_name, exon_name, exon[1], new_feature_end, feature[4]])
    new_exons_features_after.append([exon_name, exon_name, new_feature_start, exon[2], feature[4]])

    lha_end = list(filter(lambda p: 'LHA' in p, elements_list))[0][2]
    rha_start = list(filter(lambda p: 'RHA' in p, elements_list))[0][1]

    # lha_delta = lha_end - ((guide_in_transcript_pos + position_insert_start - guide_pos))
    lha_delta = lha_end - new_feature_end
    # rha_delta = rha_start - ((guide_in_transcript_pos + (len(guide) - position_insert_start - 1) - guide_pos))
    rha_delta = rha_start - new_feature_start

    for feature in new_exons_features_before:
        elements_list.append([feature[1], 
                              feature[2]+lha_delta, 
                              feature[3]+lha_delta, 
                              '+', 
                              'transcript_feature'])

    for feature in new_exons_features_after:
        elements_list.append([feature[1], 
                              feature[2]+rha_delta, 
                              feature[3]+rha_delta, 
                              '+', 
                              'transcript_feature'])
    return elements_list

title = """LOCUS       {gene_name}        {len_full_sequence} bp DNA     linear   UNA {date_today}
DEFINITION  {gene_name}.
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      synthetic DNA construct
  ORGANISM  synthetic DNA construct
REFERENCE   1  (bases 1 to {len_full_sequence})
  AUTHORS   .
  TITLE     Direct Submission
  JOURNAL   For SnapGene Viewer
            https://www.snapgene.com
FEATURES             Location/Qualifiers"""

feature_sourse = '''     source          1..{len_full_sequence}
                     /mol_type="other DNA"
                     /note="color: #ffffff"
                     /organism="synthetic DNA construct"'''

origin = '''ORIGIN
{origin_seq}
//'''

def misc_feature_template(start, end, label, color, direction):
    misc_f = f'''     misc_feature    {start}..{end}
                     /label={label}
                     /note="color: {color}; direction: {direction}"'''
    return misc_f

def primer_template(name, seq, date_today, start, end):
    primer_f = f'''     primer_bind     complement({start}..{end})
                     /label={name}
                     /note="color: black; sequence: 
                     {seq}; added: 
                     {date_today}"'''
    return primer_f

def gene_bank_file(gene_name, full_sequence, date_today, 
                   elements_list, colors, files_name, oligos = [], 
                   title=title, feature_sourse=feature_sourse, origin=origin):
    
    title = title.format(gene_name=gene_name, full_sequence=full_sequence, 
             date_today=date_today, len_full_sequence = len(full_sequence))
    
    feature_sourse = feature_sourse.format(len_full_sequence = len(full_sequence)) 
    
    all_misc_feature = ''
    for feature in elements_list:
        start = feature[1]
        end = feature[2]
        name = feature[0]
        color = colors[feature[4]][0]
        direction = feature[3]
        if direction == '+':
            direction = 'RIGHT'
        else:
            direction = 'LEFT'
        misc_feature = misc_feature_template(start, end, name, color, direction)
        all_misc_feature += misc_feature + '\n'

    all_primers = ''
    for oligo in oligos:
        end = oligo[1]
        name = oligo[0]
        seq = oligo[2]

        primer_feature = primer_template(name, seq, date_today, 1, len(seq))
        all_primers += primer_feature + '\n'
        
    origin_seq = ''
    for i in range(len(full_sequence)):
        if i%60 == 0:
            origin_seq += '\n' + ' '*(9 - len(str(i+1))) + str(i+1) + ' '
        elif i%10 == 0:
            origin_seq += ' '
        origin_seq += full_sequence[i].lower()
    origin_seq = origin_seq[1:]    
    
    origin = origin.format(origin_seq=origin_seq)

    gbk_file_name = (
        "src/static/outputs/"
        + gene_name
        + "/"
        + files_name
        + ".gbk"
    )

    with open(gbk_file_name, "w", encoding="utf-8") as f:
        f.write(title + '\n')
        f.write(feature_sourse + '\n')
        f.write(all_misc_feature)
        f.write(all_primers)
        f.write(origin + '\n')

    return gbk_file_name


# def oligo_creater(guide, full_sequence, guide_homology_arm_size, 
#                   buffer_size, n_scrambled_bases, nucleotide_changes,
#                  left_flank, right_flank, LHA_sequence, RHA_sequence,
#                  insert_sequence, elements_list, left_guide, right_guide,
#                  flank_size, is_terminal):

#     oligos = []
#     compl_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

#     left_guide_seq5, right_guide_seq5 = guide[:20], guide[-20:]

#     left_guide_homology_arm = full_sequence.split(left_guide_seq5)[0][-(guide_homology_arm_size+3):]
#     left_oligo = left_guide_homology_arm + left_guide_seq5

#     right_guide_homology_arm = full_sequence.split(right_guide_seq5)[1][:(guide_homology_arm_size+3)]
#     right_oligo = right_guide_seq5 + right_guide_homology_arm

#     if n_scrambled_bases>0:
#         left_guide_seq5_mut = (''.join([nucleotide_changes[n] for n in left_guide_seq5[:n_scrambled_bases]]) 
#                                     + left_guide_seq5[n_scrambled_bases:])
#         right_guide_seq5_mut = (right_guide_seq5[:(20-n_scrambled_bases)] 
#                                 + ''.join([nucleotide_changes[n] for n in right_guide_seq5[-n_scrambled_bases:]]))
#     else:
#         left_guide_seq5_mut = left_guide_seq5
#         right_guide_seq5_mut = right_guide_seq5

#     eff_pam = 'AGG'
#     eff_pam_rev = 'CCT'

#     left_term_oligo = right_guide_seq5 + eff_pam + LHA_sequence[:guide_homology_arm_size]
#     right_term_oligo = RHA_sequence[-guide_homology_arm_size:] + eff_pam_rev + left_guide_seq5

#     left_term_oligo = left_term_oligo.replace(right_guide_seq5, right_guide_seq5_mut)
#     right_term_oligo = right_term_oligo.replace(left_guide_seq5, left_guide_seq5_mut)
#     left_oligo = left_oligo.replace(left_guide_seq5, left_guide_seq5_mut)
#     right_oligo = right_oligo.replace(right_guide_seq5, right_guide_seq5_mut)

#     left_oligo_rev = ''.join([compl_dict[n] for n in left_oligo][::-1])
#     right_oligo_rev = ''.join([compl_dict[n] for n in right_oligo][::-1])
#     left_term_oligo_rev = ''.join([compl_dict[n] for n in left_term_oligo][::-1])
#     right_term_oligo_rev = ''.join([compl_dict[n] for n in right_term_oligo][::-1])

#     oligos.append(['LHA_adjacent_CTS', len(left_oligo_rev), left_oligo_rev])
#     oligos.append(['RHA_adjacent_CTS', len(right_oligo_rev), right_oligo_rev])

#     if is_terminal:
#         oligos.append(['LHA_terminal_CTS', len(left_term_oligo_rev), left_term_oligo_rev])
#         oligos.append(['RHA_terminal_CTS', len(right_term_oligo_rev), right_term_oligo_rev])

#         full_sequence = (left_flank[-buffer_size:] + left_term_oligo[:23] + LHA_sequence + insert_sequence 
#                             + RHA_sequence + right_term_oligo[-23:] + right_flank[:buffer_size])

#         delta  = flank_size - 23 - buffer_size
#         elements_list = [[el[i]-delta if (i==1) | (i==2) else el[i] for i in range(5)] for el in elements_list]

#     full_sequence = full_sequence.replace(right_guide_seq5, right_guide_seq5_mut)
#     full_sequence = full_sequence.replace(left_guide_seq5, left_guide_seq5_mut)

#     oligos.append(['guide', len(left_guide), left_guide])
#     oligos.append(['guide', len(right_guide), right_guide])
    
#     return full_sequence, oligos, elements_list

# def oligo_creater(guide, full_sequence, guide_homology_arm_size, 
#                   buffer_size, n_scrambled_bases, nucleotide_changes,
#                  left_flank, right_flank, LHA_sequence, RHA_sequence,
#                  insert_sequence, elements_list, left_guide, right_guide,
#                  flank_size, is_terminal):
    
def oligo_creater(guide, full_sequence, guide_homology_arm_size, 
                  buffer_size, n_scrambled_bases,
                  elements_list,
                  left_guide, right_guide, 
                  is_left_terminal, is_right_terminal, 
                  is_left, is_right):

    oligos = []
    compl_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

    # Change internal guide positions 
    atg_pos = list(filter(lambda p: p[0]=='ATG_gene', elements_list))[0]
    atg_start = atg_pos[1]-1
    atg_end = atg_pos[2]

    right_guide_start = len(full_sequence.split(right_guide)[0])
    right_guide_end = right_guide_start + len(right_guide)

    compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    left_guide_rev = ''.join([compl_dict[i] for i in left_guide][::-1])

    left_guide_start = len(full_sequence.split(left_guide_rev)[0])
    left_guide_end = left_guide_start + len(left_guide_rev)

    atg_delta_left = (left_guide_start - atg_start)%3
    if atg_delta_left==0:
        left_guide_change = full_sequence[left_guide_start-3:left_guide_end]
    elif atg_delta_left==1:
        left_guide_change = full_sequence[left_guide_start-5:left_guide_end]
    elif atg_delta_left==2:
        left_guide_change = full_sequence[left_guide_start-4:left_guide_end]

    left_guide_no_change = left_guide_change

    left_guide_change = [left_guide_change[step].upper() 
                        if step%3==2 else left_guide_change[step].lower()
                        for step in range(len(left_guide_change)) 
                        ]
    for step in range(0, len(left_guide_change)-3, 3):
        codon_i = ''.join(left_guide_change[step:step+3]).upper()
        amino_acid_i = codon_using[codon_using['CODON']==codon_i]['Amino acid'].max()
        all_codons_i = list(codon_using[(codon_using['Amino acid']==amino_acid_i)
                                        & (codon_using['CODON']!=codon_i)]['CODON'])
        if len(all_codons_i)>0:
            prefered_codon = all_codons_i[0]
            left_guide_change[step:step+3] = list(prefered_codon)

    if len(guide)==55:
        atg_delta_right = atg_delta_left
    elif len(guide)==56:
        atg_delta_right = atg_delta_left+1
    elif len(guide)==57:
        atg_delta_right = atg_delta_left+2
        if atg_delta_right==3:
            atg_delta_right = 0

    if atg_delta_right==0:
        right_guide_change = full_sequence[right_guide_start:right_guide_end+3]
    elif atg_delta_right==1:
        right_guide_change = full_sequence[right_guide_start-2:right_guide_end+3]
    elif atg_delta_right==2:
        right_guide_change = full_sequence[right_guide_start-1:right_guide_end+3]

    right_guide_no_change = right_guide_change

    right_guide_change = [right_guide_change[step].upper() 
                        if step%3==2 else right_guide_change[step].lower()
                        for step in range(len(right_guide_change)) 
                        ]
    for step in range(0, len(right_guide_change)-3, 3):
        codon_i = ''.join(right_guide_change[step:step+3]).upper()
        amino_acid_i = codon_using[codon_using['CODON']==codon_i]['Amino acid'].max()
        all_codons_i = list(codon_using[(codon_using['Amino acid']==amino_acid_i)
                                        & (codon_using['CODON']!=codon_i)]['CODON'])
        if len(all_codons_i)>0:
            prefered_codon = all_codons_i[0]
            right_guide_change[step:step+3] = list(prefered_codon)

    if is_left:
        full_sequence = full_sequence.replace(''.join(left_guide_no_change).upper(), 
                                            ''.join(left_guide_change).upper())
    if is_right:
        full_sequence = full_sequence.replace(''.join(right_guide_no_change).upper(), 
                                            ''.join(right_guide_change).upper())
    # Add additional sequences to 5' and 3' ends
    if is_left_terminal:
        full_sequence = 'AT'*4 + left_guide_no_change + full_sequence
    if is_right_terminal:
        full_sequence = full_sequence + right_guide_no_change + 'AT'*4


    oligos.append(['guide_left', len(left_guide), left_guide])
    oligos.append(['guide_right', len(right_guide), right_guide])
        
    return full_sequence, oligos

    left_guide_seq5, right_guide_seq5 = guide[:20], guide[-20:]


    left_guide_homology_arm = full_sequence.split(left_guide_seq5)[0][-(guide_homology_arm_size+3):]
    left_oligo = left_guide_homology_arm + left_guide_seq5

    right_guide_homology_arm = full_sequence.split(right_guide_seq5)[1][:(guide_homology_arm_size+3)]
    right_oligo = right_guide_seq5 + right_guide_homology_arm

    if n_scrambled_bases>0:
        left_guide_seq5_mut = (''.join([nucleotide_changes[n] for n in left_guide_seq5[:n_scrambled_bases]]) 
                                    + left_guide_seq5[n_scrambled_bases:])
        right_guide_seq5_mut = (right_guide_seq5[:(20-n_scrambled_bases)] 
                                + ''.join([nucleotide_changes[n] for n in right_guide_seq5[-n_scrambled_bases:]]))
    else:
        left_guide_seq5_mut = left_guide_seq5
        right_guide_seq5_mut = right_guide_seq5

    eff_pam = 'AGG'
    eff_pam_rev = 'CCT'

    left_term_oligo = right_guide_seq5 + eff_pam + LHA_sequence[:guide_homology_arm_size]
    right_term_oligo = RHA_sequence[-guide_homology_arm_size:] + eff_pam_rev + left_guide_seq5

    left_term_oligo = left_term_oligo.replace(right_guide_seq5, right_guide_seq5_mut)
    right_term_oligo = right_term_oligo.replace(left_guide_seq5, left_guide_seq5_mut)
    left_oligo = left_oligo.replace(left_guide_seq5, left_guide_seq5_mut)
    right_oligo = right_oligo.replace(right_guide_seq5, right_guide_seq5_mut)

    left_oligo_rev = ''.join([compl_dict[n] for n in left_oligo][::-1])
    right_oligo_rev = ''.join([compl_dict[n] for n in right_oligo][::-1])
    left_term_oligo_rev = ''.join([compl_dict[n] for n in left_term_oligo][::-1])
    right_term_oligo_rev = ''.join([compl_dict[n] for n in right_term_oligo][::-1])

    oligos.append(['LHA_adjacent_CTS', len(left_oligo_rev), left_oligo_rev])
    oligos.append(['RHA_adjacent_CTS', len(right_oligo_rev), right_oligo_rev])

    if is_terminal:
        oligos.append(['LHA_terminal_CTS', len(left_term_oligo_rev), left_term_oligo_rev])
        oligos.append(['RHA_terminal_CTS', len(right_term_oligo_rev), right_term_oligo_rev])

        full_sequence = (left_flank[-buffer_size:] + left_term_oligo[:23] + LHA_sequence + insert_sequence 
                            + RHA_sequence + right_term_oligo[-23:] + right_flank[:buffer_size])

        delta  = flank_size - 23 - buffer_size
        elements_list = [[el[i]-delta if (i==1) | (i==2) else el[i] for i in range(5)] for el in elements_list]

    full_sequence = full_sequence.replace(right_guide_seq5, right_guide_seq5_mut)
    full_sequence = full_sequence.replace(left_guide_seq5, left_guide_seq5_mut)

    oligos.append(['guide', len(left_guide), left_guide])
    oligos.append(['guide', len(right_guide), right_guide])
    
    return full_sequence, oligos, elements_list