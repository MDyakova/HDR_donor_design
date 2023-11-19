from flask import Flask, render_template, request, session, redirect, url_for, send_file
from utilities import ensemble_info, ncbi_information, guide_info, make_seqience, make_sequence_image, save_files
from flask_wtf import FlaskForm
from wtforms import StringField
import pandas as pd
import os
from io import BytesIO
import zipfile

app = Flask(__name__)
app.config['SECRET_KEY'] = 'mysecretkey'  # Change this to a secure secret key

element_sequnces_sequences = pd.read_excel('../HDR/data/all_sequences.xlsx', sheet_name='Sequences')
element_sequnces_sequences = element_sequnces_sequences.groupby(by=['Elements', 'Names'], as_index=False).max()
possible_elements = pd.unique(element_sequnces_sequences['Elements'] + '_' + element_sequnces_sequences['Names'])

out_dict = {'gene_name':'', 'ncbi_id':'', 'guide_seq':'', 'lha':400, 'rha':400, 
            'ensemble_gene_seq':'', 'gene_dict':'', 'CDS_seq':'', 'position_insert_start':'',
            'flank':100, 'guide':'', 'guide_cut_size':'', 'full_seq':'', 'possible_elements':possible_elements,
            'selected_elements':[], 'selected_elements_color':'', 'make_sequence':False,
            'all_sequences':element_sequnces_sequences, 'elements_list':[], 'full_seq_color':'',
            'image_name':'', 'insert_seq':''}

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
          'gene sequence':['#b5b5b1', '(181, 181, 177)', "<span class='grey-text'>"]}

class Gene_info(FlaskForm):
    text_field = StringField('Gene name', default='')
    text_field2 = StringField('NCBI id', default='')

class Guide_sequence(FlaskForm):
    text_field3 = StringField('Guide sequence', default='', render_kw={'style': 'width: 700px;'})

class Arm_size(FlaskForm):
    text_field4 = StringField('Size of LHA', default='', render_kw={'style': 'width: 50px;'})
    text_field5 = StringField('Size of RHA', default='', render_kw={'style': 'width: 50px;'})
    text_field6= StringField('Flank size', default='', render_kw={'style': 'width: 50px;'})



def index(out_dict, element_sequnces_sequences):

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

            out_dict['gene_name'] = gene_name
            out_dict['ncbi_id'] = ncbi_id

            ensemble_gene_seq, gene_dict = ensemble_info(gene_name)
            out_dict['ensemble_gene_seq'] = ensemble_gene_seq
            out_dict['gene_dict'] = str(gene_dict)

            transcripts_info, CDS_seq = ncbi_information(ncbi_id)
            out_dict['CDS_seq'] = CDS_seq

            # Create output directory
            output_directory = 'static/outputs/' + gene_name
            os.makedirs(output_directory, exist_ok=True)

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

            elements_list, insert_sequence, insert_sequence_color = make_seqience(out_dict['flank'], 
                                                                                out_dict['lha'] , 
                                                                                out_dict['rha'] , 
                                                                                out_dict['selected_elements'], 
                                                                                out_dict['all_sequences'],
                                                                                colors)

            out_dict['insert_seq'] = insert_sequence
            out_dict['full_seq'] = left_flank + LHA_sequence + insert_sequence + RHA_sequence + right_flank
            out_dict['full_seq_color'] = ("<span class='black-text'>" + left_flank + "</span>"
                                        + "<span class='black-text'>" + LHA_sequence + "</span>"
                                        + insert_sequence_color 
                                        + "<span class='black-text'>" + RHA_sequence + "</span>"
                                        + "<span class='black-text'>" + right_flank + "</span>")

            
            out_dict['elements_list'] = elements_list
            
            with open('out.txt', 'w') as f:
                f.write(insert_sequence_color)

        if 'del_element_submit' in request.form:
            if len(out_dict['selected_elements'])>0:
                _ = out_dict['selected_elements'].pop()
                out_dict['selected_elements_colors'] = ', '.join([colors[e.split('_')[0]][2] + e + '</span>' 
                                                                        for e in out_dict['selected_elements']])
                 
          
        selected_element = request.form.get('dropdown')
        checkbox_value = request.form.get('checkbox')
        
        if selected_element is not None:
            if checkbox_value == 'reverse':
                out_dict['selected_elements'].append(selected_element + '_reverse')
            else:
                out_dict['selected_elements'].append(selected_element)
            out_dict['selected_elements_colors'] = ', '.join([colors[e.split('_')[0]][2] + e + '</span>' 
                                                            for e in out_dict['selected_elements']])

        if 'file_upload_submit' in request.form:
            if 'file' not in request.files:
                return "No file part"
            
            file = request.files['file']
            filename = file.filename
            filename = filename.split('.')[0]
            content_type = file.content_type
            content = file.read()
            content_str = content.decode('utf-8')
            sequence = ''.join(content_str.split('\n')[1:]).strip()

            with open('out2.txt', 'w') as f:
                f.write(content_str)

            checkbox_value2 = request.form.get('checkbox2')

            if checkbox_value2 == 'in_frame':
                in_frame = 1
            else:
                in_frame = 0

            if checkbox_value == 'reverse':
                out_dict['selected_elements'].append('custom_' + filename + '_reverse')
            else:
                out_dict['selected_elements'].append('custom_' + filename)
            out_dict['selected_elements_colors'] = ', '.join([colors[e.split('_')[0]][2] + e + '</span>' 
                                                            for e in out_dict['selected_elements']])
            
            new_sequence = pd.DataFrame(data={'Elements':['custom'], 'Names':[filename], 
                   'Sequence':[sequence], 'Describe':[''],
                  'in frame': [in_frame]})
            
            out_dict['all_sequences'] = pd.concat([out_dict['all_sequences'], new_sequence])

        if 'make_image_submit' in request.form:
            make_sequence_image(out_dict['gene_name'], out_dict['elements_list'], colors, out_dict['full_seq'])
            out_dict['image_name'] = 'outputs/' + out_dict['gene_name'] + '/map_for_' + out_dict['gene_name'] + '.png'

        if 'save_files_submit' in request.form:
            fasta_file, bed_file = save_files(out_dict['gene_name'], out_dict['flank'], 
                                              out_dict['lha'], out_dict['rha'], 
                                              out_dict['elements_list'], out_dict['insert_seq'], 
                                              out_dict['full_seq'], colors)
            # return (send_file(fasta_file, as_attachment=True, download_name=fasta_file.split('/')[-1]), 
            #         send_file(bed_file, as_attachment=True, download_name=bed_file.split('/')[-1]))

            # Create a BytesIO object to store the ZIP file
            zip_buffer = BytesIO()

            # Create a ZipFile object
            with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                # Add the FASTA file to the ZIP file with a custom name
                zip_file.write(fasta_file, arcname=fasta_file.split('/')[-1])

                # Add the BED file to the ZIP file with a custom name
                zip_file.write(bed_file, arcname=bed_file.split('/')[-1])

            # Move the buffer's position to the beginning to ensure all the data is read
            zip_buffer.seek(0)

            # Return the ZIP file as an attachment
            return send_file(zip_buffer, download_name=out_dict['gene_name'] + '.zip', as_attachment=True)
           
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
    return index(out_dict, element_sequnces_sequences)

# @app.route('/button_clicked', methods=['POST'])
# def button_clicked():
#     if len(out_dict['selected_elements'])>0:
#         _ = out_dict['selected_elements'].pop()
#     return redirect(url_for('root'))

# @app.route('/button_clicked2', methods=['POST'])
# def button_clicked2():
#     out_dict['make_sequence'] = True
#     return redirect(url_for('root'))

@app.route('/download')
def download_file():

    file_path = 'readme.txt'  # Replace with the actual path to your file

    # You can customize the filename that the user will see when downloading
    filename = f'C:/Users/TargetGene/Documents/downloaded_file.txt'

    # Use send_file to send the file to the user's browser for download
    return send_file(file_path, as_attachment=True, download_name=filename)

if __name__ == '__main__':
    app.run(debug=True)