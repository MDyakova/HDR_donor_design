from flask import Flask, render_template, request, session, redirect, url_for
from utilities import ensemble_info, ncbi_information, guide_info
import pandas as pd

app = Flask(__name__)

element_sequnces_sequences = pd.read_excel('../HDR/data/all_sequences.xlsx', sheet_name='Sequences')
element_sequnces_sequences = element_sequnces_sequences.groupby(by=['Elements', 'Names'], as_index=False).max()
possible_elements = pd.unique(element_sequnces_sequences['Elements'] + '_' + element_sequnces_sequences['Names'])

out_dict = {'gene_name':'', 'ncbi_id':'', 'guide_seq':'', 'lha':'', 'rha':'', 
            'ensemble_gene_seq':'', 'gene_dict':'', 'CDS_seq':'', 'position_insert_start':'',
            'flank':'', 'guide':'', 'guide_cut_size':'', 'full_seq':'', 'possible_elements':possible_elements,
            'selected_elements':[]}



def index(out_dict):
    
    if request.method == 'POST':
        gene_name = request.form.get('gene_name')
        ncbi_id = request.form.get('ncbi_id')
        guide_seq = request.form.get('text1')

        lha = request.form.get('lha')
        rha = request.form.get('rha')
        flank = request.form.get('flank')

        if gene_name is not None:
            out_dict['gene_name'] = gene_name
        if ncbi_id is not None:
            out_dict['ncbi_id'] = ncbi_id
        if guide_seq is not None:
            out_dict['guide_seq'] = guide_seq
        if lha is not None:
            out_dict['lha'] = int(lha)
        if rha is not None:
            out_dict['rha'] = int(rha)
        if flank is not None:
            out_dict['flank'] = int(flank)

        if gene_name is not None:
            ensemble_gene_seq, gene_dict = ensemble_info(gene_name)
            out_dict['ensemble_gene_seq'] = ensemble_gene_seq
            out_dict['gene_dict'] = str(gene_dict)

        if ncbi_id is not None:
            transcripts_info, CDS_seq = ncbi_information(ncbi_id)
            out_dict['CDS_seq'] = CDS_seq

        if guide_seq is not None:
            position_insert_start, guide_cut_size, guide = guide_info(guide_seq, out_dict['CDS_seq'])
            out_dict['position_insert_start'] = position_insert_start
            out_dict['guide'] = guide
            out_dict['guide_cut_size'] = guide_cut_size

        if out_dict['lha']!='':
            LHA_sequence = (out_dict['ensemble_gene_seq'].split(out_dict['guide'])[0][-(out_dict['lha'] 
                                                                           - out_dict['position_insert_start']):] 
                                                                           + out_dict['guide'][:out_dict['position_insert_start']])
            RHA_sequence = (out_dict['guide'][-out_dict['guide_cut_size']:] 
                            + out_dict['ensemble_gene_seq'].split(out_dict['guide'])[1][:(out_dict['rha'] - out_dict['guide_cut_size'])])
            left_flank = out_dict['ensemble_gene_seq'].split(LHA_sequence)[0][-out_dict['flank']:]
            right_flank = out_dict['ensemble_gene_seq'].split(RHA_sequence)[1][:out_dict['flank']]

            insert_sequence = ''

            # out_dict['full_seq'] = left_flank + LHA_sequence + insert_sequence + RHA_sequence + right_flank
            out_dict['full_seq'] = ("<span class='black-text'>" + left_flank 
                                    + "</span><span class='violet-text'>" + LHA_sequence 
                                    + "</span><span class='black-text'>" + insert_sequence 
                                    + "</span><span class='violet-text'>" + RHA_sequence 
                                    + "</span><span class='black-text'>" + right_flank 
                                    + "</span>.")
            
        selected_element = request.form.get('dropdown')
        if selected_element is not None:
            out_dict['selected_elements'].append(selected_element)
            print(out_dict['selected_elements'])

        item_to_delete = request.form.get('delete_item')
        print(item_to_delete)


        return render_template('home.html', **out_dict)

    return render_template('home.html', **out_dict)

@app.route('/', methods=['GET', 'POST'])
def root():
    return index(out_dict)

@app.route('/button_clicked', methods=['POST'])
def button_clicked():
    if len(out_dict['selected_elements'])>0:
        _ = out_dict['selected_elements'].pop()
    return redirect(url_for('root'))

if __name__ == '__main__':
    app.run(debug=True)