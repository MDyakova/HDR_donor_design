import pandas as pd
import numpy as np
from pyensembl import EnsemblRelease
import matplotlib.pyplot as plt
from dna_features_viewer import GraphicFeature, GraphicRecord
from Bio import Entrez
from Bio import SeqIO
from matplotlib import colors as m_colors
from datetime import date
import requests
from matplotlib import colors as m_colors

def ensemble_info(gene_name):
    # release 109 uses human reference genome GRCh38
    data = EnsemblRelease(109)

    data.download() 
    data.index()

    gene = data.genes_by_name(gene_name)[0]
    exon_ids  = data.exon_ids_of_gene_name(gene_name)
    transcripts_id = data.transcript_ids_of_gene_name(gene_name)
    strand = gene.strand
    chromosome_name = gene.contig

    gene_dict = gene.to_dict()

    # if os.path.exists('outputs/' + gene_name) is not True:
    #     os.mkdir('outputs/' + gene_name)
        
    # # https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/
    # with open('../gene_information/databases/reference_genome/Homo_sapiens.GRCh38.dna.chromosome.' + chromosome_name + '.fa') as f:
    #     seq = f.readlines()
        
    # seq = [i.replace('\n', '') for i in seq[1:]]
    # full_sequence = ''.join(seq)


    start_gene = gene.start - 1500
    end_gene = gene.end + 1500

    url = f'https://rest.ensembl.org/sequence/region/human/{chromosome_name}:{start_gene}-{end_gene}?content-type=application/json;version=109'
    data = requests.get(url)
    gene_seq = data.json()['seq']

    # gene_seq = list(full_sequence[(start_gene - 1):(end_gene)])
    coord_list =  list([i for i in range(start_gene, end_gene + 1)])

    if strand == '-':
        compl_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
        gene_seq = ([compl_dict[l] for l in gene_seq][::-1])
        coord_list = coord_list[::-1]
        
    ensemble_gene_seq = ''.join(gene_seq)

    del gene_dict['gene_name']
    del gene_dict['gene_id']
    del gene_dict['genome']

    return ensemble_gene_seq, gene_dict

def ncbi_information(NCBI_id):
    # RefSeq information

    transcripts_info = []
    Entrez.email = "Your.Name.Here@example.org"

    handle = Entrez.efetch(db="nucleotide", id=NCBI_id, rettype="gb", retmode="text")
    seq_record = [seq_record for seq_record in SeqIO.parse(handle, "gb")][0]
    refseq_sequence = str(seq_record.seq)
    features_list = []
    for feature in seq_record.features:
        features_list.append([feature.type, int(feature.location.start+1), int(feature.location.end), 
                            feature.qualifiers])
        if feature.type == 'CDS':
            protein_id = feature.qualifiers['protein_id'][0]

    annotation = seq_record.annotations['keywords']
    ensemble_is_rs = ''
    if 'MANE Select' in annotation:
        ensemble_is_rs = seq_record.annotations['structured_comment']['RefSeq-Attributes']['MANE Ensembl match'].split('/')[0].split('.')[0]

    transcripts_info.append([NCBI_id, refseq_sequence, features_list, ensemble_is_rs, protein_id])    
    handle.close()

    for rna_info in transcripts_info:
        refseq_s = rna_info[1]
        features = rna_info[2]
        ensemble_is_rs = rna_info[3]
        protein_id = rna_info[4]
        CDS = list(filter(lambda p: p[0] == 'CDS', features))[0]
        CDS_start_pos = CDS[1]
        CDS_end_pos = CDS[2]

    CDS_seq = refseq_sequence[CDS_start_pos-1:CDS_end_pos]

    return transcripts_info, CDS_seq

def guide_info(guide_seq, CDS_seq):
    compl_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    guide = ''.join([compl_dict[i] for i in guide_seq][::-1])

    guide_cut_size = len(guide)//2

    if guide[:guide_cut_size] in CDS_seq:
        guide_in_CDS_seq = guide[:guide_cut_size]
        guide_in_CDS_pos = len(CDS_seq.split(guide[:guide_cut_size])[0])
        guide_pos = 0
    elif guide[-guide_cut_size:] in CDS_seq:
        guide_in_CDS_seq = guide[-guide_cut_size:]
        guide_in_CDS_pos = len(CDS_seq.split(guide[-guide_cut_size:])[0])
        guide_pos = len(guide.split(guide[-guide_cut_size:])[0])
    else:
        print('guide not found')

    codon_positions = {}
    for i, k in zip(range(len(guide)), range(len(guide))):
        codon_positions[k] = (guide_in_CDS_pos - (guide_pos - i))%3
        
    cut_site_codon_pos = codon_positions[guide_cut_size]
    position_insert_start = guide_cut_size - cut_site_codon_pos

    return position_insert_start, guide_cut_size, guide


def make_seqience(flank_size, LHA_size, RHA_size, donor_elements, element_sequnces_sequences, colors):

    stop_codons = ['TAA', 'TAG', 'TGA']
    compl_dict = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}

    insert_start = flank_size + LHA_size + 1
    next_element = insert_start

    elements_list = []
    elements_list.append(['LHA', flank_size+1, flank_size + LHA_size, '+', 'gene sequence'])

    insert_sequence = ''
    insert_sequence_color = ''
    coding_sequence = ''
    new_start_codon = False
    for step, element in enumerate(donor_elements):
        element = '_'.join(element.split('_')[1:])
        if '_reverse' in element:
            direction = '-'
            element = element.split('_reverse')[0]
            seq_i = element_sequnces_sequences[element_sequnces_sequences['Names']==element]['Sequence'].max().upper()
            group = element_sequnces_sequences[element_sequnces_sequences['Names']==element]['Elements'].max()
            seq_i = seq_i.replace('U', 'T')
            seq_i = ''.join([compl_dict[i] for i in seq_i][::-1])
        else:
            direction = '+'
            seq_i = element_sequnces_sequences[element_sequnces_sequences['Names']==element]['Sequence'].max().upper()
            group = element_sequnces_sequences[element_sequnces_sequences['Names']==element]['Elements'].max()
            seq_i = seq_i.replace('U', 'T')
        
            
        in_frame = element_sequnces_sequences[element_sequnces_sequences['Names']==element]['in frame'].max()
        if (seq_i[-3:] in stop_codons) & (in_frame==1):
            seq_i = seq_i[:-3]
            print(element, 'STOP')
        
        insert_sequence += seq_i

        if elements_list[-1][4] == 'Promoter':
            coding_sequence = ''
            new_start_codon = True
        coding_sequence += seq_i

        insert_sequence_color += colors[group][2]
        insert_sequence_color += seq_i
        insert_sequence_color += "</span>"

        # if (len(insert_sequence)%3 != 0) & (in_frame==1):
        #     insert_sequence += 'N'*(3 - len(insert_sequence)%3)
        #     insert_sequence_color += colors[group][2]
        #     insert_sequence_color += 'N'*(3 - len(insert_sequence)%3)
        #     insert_sequence_color += "</span>"

        if (len(coding_sequence)%3 != 0) & (in_frame==1):
            if new_start_codon:
                atg_ind = np.argmax([coding_sequence[i:i+3]=='ATG' for i in range(len(coding_sequence))])
                new_start_codon = False
                if len(coding_sequence[atg_ind:])%3!=0:
                    insert_sequence += 'N'*(3 - len(coding_sequence[atg_ind:])%3)
                    coding_sequence += 'N'*(3 - len(coding_sequence[atg_ind:])%3)
                    insert_sequence_color += colors[group][2]
                    insert_sequence_color += 'N'*(3 - len(coding_sequence[atg_ind:])%3)
                    insert_sequence_color += "</span>"
            else:
                insert_sequence += 'N'*(3 - len(insert_sequence)%3)
                insert_sequence_color += colors[group][2]
                insert_sequence_color += 'N'*(3 - len(insert_sequence)%3)
                insert_sequence_color += "</span>"

           
        # print(seq_i)
        # print()
        # print(insert_sequence_color)

        elements_list.append([element, next_element, insert_start + len(insert_sequence) - 1, direction, group])
        next_element = insert_start + len(insert_sequence)

    elements_list.append(['RHA', elements_list[-1][2]+1, elements_list[-1][2] + RHA_size, '+', 'gene sequence'])

    return elements_list, insert_sequence, insert_sequence_color

def make_sequence_image(gene_name, elements_list, colors, full_sequence):
    features = []
    plt.figure(figsize=(10, 20))
    for element in elements_list:
        if element[3] == '+':
            strand_plot = +1
        else:
            strand_plot = -1
        feature = GraphicFeature(start=element[1], end=element[2], strand=strand_plot, color=colors[element[4]][0],
                            label=element[0])
        features.append(feature)
    #     break
    record = GraphicRecord(sequence_length=len(full_sequence), features=features)
    cropped_record = record.crop((1, len(full_sequence)))
    cropped_record.sequence = full_sequence
    cropped_record.plot(figure_width=10, strand_in_label_threshold=7, plot_sequence=False)
    plt.savefig('src/static/outputs/' + gene_name + '/map_for_' + gene_name + '.png', bbox_inches='tight')


def save_files(gene_name, flank_size, LHA_size, RHA_size, elements_list, insert_sequence, full_sequence, colors):
    date_today = str(date.today())
    features_list_donor = []

    # features_list_donor.append(['0', int(flank_size + 1), int(flank_size + LHA_size), 'LHA', 1000, '+',
    #                                 0, 1, ','.join([str(int(c*255)) for c in m_colors.to_rgb(colors['gene sequence'][0])]), 
    #                                 None, None, None])
    # features_list_donor.append(['0', int(flank_size + LHA_size + len(insert_sequence) + 1), 
    #                                 int(flank_size + LHA_size + len(insert_sequence) + RHA_size), 'RHA', 1000, '+',
    #                                 0, 1, ','.join([str(int(c*255)) for c in m_colors.to_rgb(colors['gene sequence'][0])]), 
    #                                 None, None, None])

    print(elements_list)
    for feature in elements_list:
        start = feature[1]
        end = feature[2]
        name = feature[0]
        color = colors[feature[4]][0]
        direction = feature[3]
        features_list_donor.append(['0', int(start), int(end), name, 1000, direction,
                                    0, 1, ','.join([str(int(c*255)) for c in m_colors.to_rgb(color)]), None, None, None])
        
    features_df = pd.DataFrame(features_list_donor)
    features_df[1] = features_df[1] - 1

    features_df.to_csv('static/outputs/' + gene_name + '/' + gene_name + '_donor.bed', sep='\t', header=None, index=None)

    with open('static/outputs/' + gene_name + '/' + gene_name + '_donor.bed', 'r') as f:
        bed_file = f.read()

    bed_file_name = 'static/outputs/' + gene_name + '/' + gene_name + '_donor_' + date_today + '.bed'
    with open(bed_file_name, 'w') as f:
        f.write('track itemRgb=On\n')
        f.write(bed_file)

    fasta_file_name = 'static/outputs/' + gene_name + '/' + gene_name + '_donor_sequence.fa'
    with open(fasta_file_name, 'w') as f:
        f.write('> ' + gene_name + '_donor_sequence' + '\n')
        f.write(full_sequence + '\n')

    return fasta_file_name, bed_file_name