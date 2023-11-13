from flask import Flask, render_template, request, session, redirect, url_for
from utilities import ensemble_info, ncbi_information, guide_info, make_seqience
from flask_wtf import FlaskForm
from wtforms import StringField
import pandas as pd

app = Flask(__name__)
app.config['SECRET_KEY'] = 'mysecretkey'  # Change this to a secure secret key

element_sequnces_sequences = pd.read_excel('../HDR/data/all_sequences.xlsx', sheet_name='Sequences')
element_sequnces_sequences = element_sequnces_sequences.groupby(by=['Elements', 'Names'], as_index=False).max()
possible_elements = pd.unique(element_sequnces_sequences['Elements'] + '_' + element_sequnces_sequences['Names'])

out_dict = {'gene_name':'', 'ncbi_id':'', 'guide_seq':'', 'lha':400, 'rha':400, 
            'ensemble_gene_seq':'', 'gene_dict':'', 'CDS_seq':'', 'position_insert_start':'',
            'flank':100, 'guide':'', 'guide_cut_size':'', 'full_seq':'', 'possible_elements':possible_elements,
            'selected_elements':[], 'make_sequence':False}

colors = {'2A motif':'#0000EE',
          'protein': '#00C957',
          'cloning' : '#CDB38B',
          'Stop codon': '#FF3030',
          'Terminator': '#FF6103'}

class Gene_info(FlaskForm):
    text_field = StringField('Gene name', default='')
    text_field2 = StringField('NCBI id', default='')

class Guide_sequence(FlaskForm):
    text_field3 = StringField('Guide sequence', default='', render_kw={'style': 'width: 700px;'})

class Arm_size(FlaskForm):
    text_field4 = StringField('Size of LHA', default='', render_kw={'style': 'width: 50px;'})
    text_field5 = StringField('Size of RHA', default='', render_kw={'style': 'width: 50px;'})
    text_field6= StringField('Flank size', default='', render_kw={'style': 'width: 50px;'})



def index(out_dict):

    gene_info_form = Gene_info()
    guide_seq_form = Guide_sequence()
    arm_size_form = Arm_size()

    arm_size_form.text_field4.default = out_dict['lha']
    arm_size_form.text_field5.default = out_dict['rha']
    arm_size_form.text_field6.default = out_dict['flank']
    arm_size_form.process()

    forms={'gene_info_form':gene_info_form, 
           'guide_seq_form':guide_seq_form, 
           'arm_size_form':arm_size_form}


    if request.method == 'POST':

        if 'gene_info_form_submit' in request.form:
            # Process form submission if needed
            gene_name = gene_info_form.text_field.data
            ncbi_id = gene_info_form.text_field2.data

            print(ncbi_id)

            out_dict['gene_name'] = gene_name
            out_dict['ncbi_id'] = ncbi_id

            ensemble_gene_seq, gene_dict = ensemble_info(gene_name)
            out_dict['ensemble_gene_seq'] = ensemble_gene_seq
            out_dict['gene_dict'] = str(gene_dict)

            transcripts_info, CDS_seq = ncbi_information(ncbi_id)
            out_dict['CDS_seq'] = CDS_seq

        if 'guide_seq_form_submit' in request.form:
            guide_seq = guide_seq_form.text_field3.data


            out_dict['guide_seq'] = guide_seq

            position_insert_start, guide_cut_size, guide = guide_info(guide_seq, out_dict['CDS_seq'])
            out_dict['position_insert_start'] = position_insert_start
            out_dict['guide'] = guide
            out_dict['guide_cut_size'] = guide_cut_size

        if 'make_seq_submit' in request.form:
            out_dict['make_sequence'] = True

            lha = arm_size_form.text_field4.data
            rha = arm_size_form.text_field5.data
            flank = arm_size_form.text_field6.data

            out_dict['lha'] = int(lha)
            out_dict['rha'] = int(rha)
            out_dict['flank'] = int(flank)
    
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
            
        if 'del_element_submit' in request.form:
            if len(out_dict['selected_elements'])>0:
                _ = out_dict['selected_elements'].pop()

    
          
        selected_element = request.form.get('dropdown')
        checkbox_value = request.form.get('checkbox')
        
        if selected_element is not None:
            if checkbox_value == 'reverse':
                out_dict['selected_elements'].append(selected_element + '_reverse')
            else:
                out_dict['selected_elements'].append(selected_element)

        item_to_delete = request.form.get('delete_item')

        
        gene_info_form.text_field.default = out_dict['gene_name']
        gene_info_form.text_field2.default = out_dict['ncbi_id']
        gene_info_form.process()

        guide_seq_form.text_field3.default = out_dict['guide_seq']
        guide_seq_form.process()
        
        
        forms={'gene_info_form':gene_info_form, 
                'guide_seq_form':guide_seq_form, 
                'arm_size_form':arm_size_form}

        return render_template('home.html', out_dict = out_dict, forms=forms)

    return render_template('home.html', out_dict = out_dict, forms=forms)

@app.route('/', methods=['GET', 'POST'])
def root():
    return index(out_dict)

# @app.route('/button_clicked', methods=['POST'])
# def button_clicked():
#     if len(out_dict['selected_elements'])>0:
#         _ = out_dict['selected_elements'].pop()
#     return redirect(url_for('root'))

# @app.route('/button_clicked2', methods=['POST'])
# def button_clicked2():
#     out_dict['make_sequence'] = True
#     return redirect(url_for('root'))

if __name__ == '__main__':
    app.run(debug=True)