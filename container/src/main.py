"""
Code launch FLASK server to make HDR donor template
"""
import os
from io import BytesIO
import zipfile
import pandas as pd
import numpy as np
from datetime import date
import json

from flask import Flask, render_template, request, send_file
from flask_wtf import FlaskForm
from wtforms import StringField

from utilities import (
    ensemble_info,
    ncbi_information,
    guide_info,
    make_seqience,
    make_sequence_image,
    save_files,
    list_files_and_sizes,
    find_promoter, 
    gene_features,
    gene_bank_file,
    oligo_creater
)

app = Flask(__name__)
app.config["SECRET_KEY"] = "mysecretkey"  # fake key to work with flask server

# table with all possible elements and sequences
element_sequences = pd.read_excel(
    "src/static/data/all_sequences.xlsx", sheet_name="Sequences"
)
element_sequences = element_sequences.groupby(
    by=["Elements", "Names"], as_index=False
).max()
element_sequences = element_sequences.sort_values(
    by=["Elements", "Names"], key=lambda col: col.str.lower()
)
possible_elements = [
    {"group": i, "name": j}
    for i, j in zip(
        element_sequences["Elements"], element_sequences["Names"]
    )
]

# dictionary with all useful variables necessary to save throw whole pipeline
with open('config.json', 'r') as f:
    config = json.load(f)

out_dict = config['initial_values']
out_dict['possible_elements'] = possible_elements
out_dict['all_sequences'] = element_sequences
out_dict['make_sequence'] = False

out_dict_start = out_dict.copy()

# colors for different elements
colors = config['colors']

class GeneInfo(FlaskForm):
    """
    Information about ensemble and ncbi gene names or id
    """
    text_field = StringField("Ensemble gene name", default="")
    text_field2 = StringField("  NCBI id", default="")
    text_field3 = StringField(
        "Guide sequence", default="", render_kw={"style": "width: 550px;"}
    )
    text_field4 = StringField(
        "Size of LHA", default="", render_kw={"style": "width: 50px;"}
    )
    text_field5 = StringField(
        "Size of RHA", default="", render_kw={"style": "width: 50px;"}
    )
    text_field6 = StringField(
        "Primers place size", default="", render_kw={"style": "width: 50px;"}
    )

class SaveFiles(FlaskForm):
    """
    Information about ensemble and ncbi gene names or id
    """
    text_field7 = StringField("", default="", render_kw={"style": "width: 550px;"})

class CTSInfo(FlaskForm):
    """
    Information about ensemble and ncbi gene names or id
    """
    text_field8 = StringField("Size of CTS homology arm", default="", render_kw={"style": "width: 50px;"})
    text_field9 = StringField("Size of buffer", default="", render_kw={"style": "width: 50px;"})
    text_field10 = StringField("Number of scrambled bases", default="", render_kw={"style": "width: 50px;"})



def index(out_dict):
    """
    Launch main code to prepare input data to donor sequence
    """

    gene_info_form = GeneInfo()
    save_files_form = SaveFiles()
    cts_info_form = CTSInfo()

    forms = {
        "gene_info_form": gene_info_form,
        "save_files_form":save_files_form,
        "cts_info_form":cts_info_form

    }

    if request.method == "POST":
        if "gene_info_form_submit" in request.form:
            # Data from ensemble and ncbi, guide sequence
            gene_name = gene_info_form.text_field.data
            gene_name = gene_name.upper()
            ncbi_id = gene_info_form.text_field2.data
            ncbi_id = ncbi_id.upper()
            guide_seq = gene_info_form.text_field3.data
            guide_seq = guide_seq.upper()

            if (gene_name == '') | (ncbi_id == '') | (guide_seq == ''):
                text_error = 'enter all data'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms)


            out_dict["gene_name"] = gene_name
            out_dict["ncbi_id"] = ncbi_id

            try:
                ensemble_gene_seq, gene_dict, strand = ensemble_info(gene_name)
                out_dict["ensemble_gene_seq"] = ensemble_gene_seq
                out_dict["strand"] = strand
                out_dict["gene_dict"] = 'Gene info: ' + str(gene_dict)
            except Exception as e:
                text_error = 'check gene name'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms)

            try:
                _, cds_seq, cds_seq_long, features_list, refseq_sequence = ncbi_information(ncbi_id)
                out_dict["CDS_seq"] = cds_seq
                out_dict["CDS_seq_long"] = cds_seq_long
                out_dict["features_list"] = features_list
                out_dict["refseq_seq"] = refseq_sequence
            except Exception as e:
                text_error = 'check NCBI id'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms)

            # Create output directory
            output_directory = "src/static/outputs/" + gene_name
            os.makedirs(output_directory, exist_ok=True)

            # guide sequence and position
            # guide_seq = gene_info_form.text_field3.data

            out_dict["guide_seq"] = guide_seq

            try:
                (position_insert_start, 
                 guide_cut_size, 
                 guide, 
                 guide_in_cds_pos, 
                 guide_in_cds_seq, 
                 guide_in_transcript_pos, 
                 guide_pos, left_guide, right_guide) = guide_info(
                    guide_seq, out_dict["CDS_seq"], out_dict["strand"], 
                    out_dict["ensemble_gene_seq"], out_dict["CDS_seq_long"], out_dict["refseq_seq"]
                )

                out_dict["position_insert_start"] = position_insert_start
                out_dict["guide"] = guide
                out_dict["guide_cut_size"] = guide_cut_size
                out_dict["guide_in_cds_pos"] = guide_in_cds_pos
                out_dict["guide_in_cds_seq"] = guide_in_cds_seq
                out_dict["guide_in_transcript_pos"] = guide_in_transcript_pos
                out_dict["guide_pos"] = guide_pos
                out_dict["left_guide"] = left_guide
                out_dict["right_guide"] = right_guide

            except Exception as e:
                text_error = 'check guide sequence'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms)

            # arm and flank sizes
            lha = gene_info_form.text_field4.data
            rha = gene_info_form.text_field5.data
            flank = gene_info_form.text_field6.data

            try:
                out_dict["lha"] = int(lha)
                out_dict["rha"] = int(rha)
                out_dict["flank"] = int(flank)
            except Exception as e:
                text_error = 'arms and flank sizes must be > 0'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms)
            date_today = str(date.today())
            file_names = ('Donor_' 
                            + out_dict["gene_name"] 
                            + '_' 
                            + out_dict["guide_seq"]
                            + '_'
                            + date_today)
            out_dict['files_name'] = file_names

        if "make_seq_submit" in request.form:
            # data to make full donor sequence
            out_dict["make_sequence"] = True

            lha_sequence = (
                out_dict["ensemble_gene_seq"].split(out_dict["guide"])[0][
                    -(out_dict["lha"] - out_dict["position_insert_start"]) :
                ]
                + out_dict["guide"][: out_dict["position_insert_start"]]
            )
            rha_sequence = (
                out_dict["guide"][-out_dict["guide_cut_size"] :]
                + out_dict["ensemble_gene_seq"].split(out_dict["guide"])[1][
                    : (out_dict["rha"] - out_dict["guide_cut_size"])
                ]
            )

            delta_nucleotides = len(out_dict["guide"].split(out_dict["guide"][:out_dict["position_insert_start"]])[1].split(out_dict["guide"][-out_dict["guide_cut_size"]:])[0])
            out_dict['delta_nucleotides'] = delta_nucleotides

            left_flank = out_dict["ensemble_gene_seq"].split(lha_sequence)[0][
                -out_dict["flank"] :
            ]
            right_flank = out_dict["ensemble_gene_seq"].split(rha_sequence)[1][
                : out_dict["flank"]
            ]

            out_dict['left_flank'] = left_flank
            out_dict['right_flank'] = right_flank
            out_dict['lha_sequence'] = lha_sequence
            out_dict['rha_sequence'] = rha_sequence

            promoter_list = find_promoter(out_dict["gene_name"])
            # left_seq = left_flank + lha_sequence + rha_sequence[:20]
            left_seq = (left_flank 
                        + out_dict["ensemble_gene_seq"].split(out_dict["guide"])[0][-(out_dict["lha"] - out_dict["position_insert_start"]):] 
                        + out_dict["guide"])

            insert_sequence = ""

            elements_list, insert_sequence, insert_sequence_color = make_seqience(
                out_dict["flank"],
                out_dict["lha"],
                out_dict["rha"],
                out_dict["selected_elements"],
                out_dict["all_sequences"],
                colors,
                out_dict["guide_in_cds_pos"],
                promoter_list,
                left_seq,
                out_dict["CDS_seq"]

            )

            if 'N' in insert_sequence:
                out_dict["nn_error"] = ("<span class='red-text'>" 
                                         + 'Attention! Frameshift in coding sequence. Added "N" nucleotides.'
                                         + "</span>")

            out_dict["insert_seq"] = insert_sequence
            out_dict["full_seq"] = (
                left_flank + lha_sequence + insert_sequence + rha_sequence + right_flank
            )
            out_dict["full_seq_color"] = (
                "<span class='black-text'>"
                + left_flank
                + "</span>"
                + "<span class='black-text'>"
                + lha_sequence
                + "</span>"
                + insert_sequence_color
                + "<span class='black-text'>"
                + rha_sequence
                + "</span>"
                + "<span class='black-text'>"
                + right_flank
                + "</span>"
            )

            if np.max(['exon' in i[0] for i in out_dict["features_list"]]):
                elements_list = gene_features(out_dict["features_list"], 
                                            out_dict["guide"], 
                                            out_dict["guide_pos"], 
                                            out_dict["guide_in_transcript_pos"], 
                                            elements_list,
                                            out_dict["position_insert_start"],
                                            out_dict['delta_nucleotides'])
            
            out_dict["elements_list"] = elements_list

        if "del_element_submit" in request.form:
            # delete selected element
            if len(out_dict["selected_elements"]) > 0:
                _ = out_dict["selected_elements"].pop()
                out_dict["selected_elements_colors"] = ", ".join(
                    [
                        colors[e.split("_")[0]][2] + e + "</span>"
                        for e in out_dict["selected_elements"]
                    ]
                )

        selected_element = request.form.get("dropdown")
        checkbox_value = request.form.get("checkbox")

        if selected_element is not None:
            if selected_element != '':
                if checkbox_value == "reverse":
                    out_dict["selected_elements"].append(selected_element + "_reverse")
                else:
                    out_dict["selected_elements"].append(selected_element)
                out_dict["selected_elements_colors"] = ", ".join(
                    [
                        colors[e.split("_")[0]][2] + e + "</span>"
                        for e in out_dict["selected_elements"]
                    ]
                )

        if "file_upload_submit" in request.form:
            # upload fasta file with new element
            if "file" not in request.files:
                return "No file part"

            file = request.files["file"]
            filename = file.filename
            filename = filename.split(".")[0]
            # content_type = file.content_type
            content = file.read()
            content_str = content.decode("utf-8")
            sequence = "".join(content_str.split("\n")[1:]).replace('\r', '').strip()

            # with open('out2.txt', 'w') as f:
            #     f.write(content_str)

            checkbox_value2 = request.form.get("checkbox2")
            checkbox_value3 = request.form.get("checkbox3")

            if checkbox_value2 == "in_frame":
                in_frame = 1
            else:
                in_frame = 0

            if checkbox_value3 == "reverse":
                out_dict["selected_elements"].append("custom_" + filename + "_reverse")
            else:
                out_dict["selected_elements"].append("custom_" + filename)
            out_dict["selected_elements_colors"] = ", ".join(
                [
                    colors[e.split("_")[0]][2] + e + "</span>"
                    for e in out_dict["selected_elements"]
                ]
            )

            new_sequence = pd.DataFrame(
                data={
                    "Elements": ["custom"],
                    "Names": [filename],
                    "Sequence": [sequence],
                    "Describe": [""],
                    "in frame": [in_frame],
                }
            )

            out_dict["all_sequences"] = pd.concat(
                [out_dict["all_sequences"], new_sequence]
            )

        if "make_image_submit" in request.form:
            # make image with all selected elements
            make_sequence_image(
                out_dict["gene_name"],
                out_dict["elements_list"],
                colors,
                out_dict["full_seq"],
            )
            out_dict["image_name"] = (
                "outputs/"
                + out_dict["gene_name"]
                + "/map_for_"
                + out_dict["gene_name"]
                + ".png"
            )

        if "save_files_submit" in request.form:
            # save files for SnapGene
            # files_name = request.form['save_files_input']
            files_name = save_files_form.text_field7.data
            out_dict["files_name"] = files_name

            fasta_file, bed_file = save_files(
                out_dict["gene_name"],
                out_dict["elements_list"],
                out_dict["full_seq"],
                colors,
                out_dict["files_name"]
            )

            date_today = str(date.today())
            gbk_file = gene_bank_file(out_dict["gene_name"], out_dict["full_seq"], date_today, 
                                            out_dict["elements_list"], colors, out_dict["files_name"])

            # Create a BytesIO object to store the ZIP file
            zip_buffer = BytesIO()

            # Create a ZipFile object
            with zipfile.ZipFile(
                zip_buffer, "a", zipfile.ZIP_DEFLATED, False
            ) as zip_file:
                # Add the FASTA file to the ZIP file with a custom name
                zip_file.write(fasta_file, arcname=fasta_file.split("/")[-1])

                # Add the BED file to the ZIP file with a custom name
                zip_file.write(bed_file, arcname=bed_file.split("/")[-1])

                # Add the GBK file to the ZIP file with a custom name
                zip_file.write(gbk_file, arcname=gbk_file.split("/")[-1])

            # Move the buffer's position to the beginning to ensure all the data is read
            zip_buffer.seek(0)

            # Return the ZIP file as an attachment
            return send_file(
                zip_buffer,
                download_name=out_dict["files_name"] + ".zip",
                as_attachment=True,
            )
        
        if "clear_forms_submit" in request.form:
            for key in out_dict.keys():
                out_dict[key] = out_dict_start[key]
                out_dict["selected_elements"] = []

        if "cts_info_form_submit" in request.form:
            # Data from ensemble and ncbi, guide sequence
            cts_ha_size = int(cts_info_form.text_field8.data)
            buffer = int(cts_info_form.text_field9.data)
            scrambled_nt = int(cts_info_form.text_field10.data)

            if (cts_ha_size == '') | (buffer == '') | (scrambled_nt == ''):
                text_error = 'enter all data'
                out_dict["gene_dict"] = ("<span class='red-text'>" 
                                         + 'Error: ' + str(text_error)
                                         + "</span>")
                return render_template("home.html", out_dict=out_dict, forms=forms)

            out_dict["CTS_HA"] = np.maximum(cts_ha_size, 23)
            out_dict["buffer"] = buffer
            out_dict["scrambled_nt"] = scrambled_nt           

            # files_name = save_files_form.text_field7.default
            # out_dict["files_name"] = files_name

            nucleotide_changes = {'A':'G', 'T':'C', 'C':'T', 'G':'A'}

            checkbox_value4 = request.form.get("checkbox4")

            if checkbox_value4 == "is_terminal_oligos":
                is_terminal = True
            else:
                is_terminal = False

            full_sequence, oligos, elements_list = oligo_creater(out_dict["guide"], out_dict["full_seq"], out_dict["CTS_HA"], 
                                                                out_dict["buffer"], out_dict["scrambled_nt"], nucleotide_changes,
                                                                out_dict['left_flank'], out_dict['right_flank'], 
                                                                out_dict['lha_sequence'], out_dict['rha_sequence'],
                                                                out_dict["insert_seq"], out_dict["elements_list"],
                                                                out_dict['left_guide'], out_dict['right_guide'],
                                                                out_dict["flank"], is_terminal)
            
            out_dict["elements_list_oligo"] = elements_list
            out_dict["full_seq_oligo"] = full_sequence
            out_dict["oligos"] = oligos

            date_today = str(date.today())
            gbk_file = gene_bank_file(out_dict["gene_name"], out_dict["full_seq_oligo"], date_today, 
                                            out_dict["elements_list_oligo"], colors, out_dict["files_name"], 
                                            oligos = out_dict["oligos"])
            
            fasta_file_name = (
                "src/static/outputs/" + out_dict["gene_name"] + "/" + out_dict["files_name"] + "_donor_sequence_with_oligo.fa"
            )
            with open(fasta_file_name, "w", encoding="utf-8") as file:
                file.write("> " + out_dict["gene_name"] + "\n")
                file.write(out_dict["full_seq_oligo"] + "\n")

            fasta_file_name_oligos = (
                            "src/static/outputs/" + out_dict["gene_name"] + "/" + out_dict["files_name"] + "_oligo_sequences.fa"
                        )

            with open(fasta_file_name_oligos, "w", encoding="utf-8") as file:
                for oligo in out_dict["oligos"]:
                    if 'guide' not in oligo[0]:
                        file.write('> ' + oligo[0] + '\n')
                        file.write(oligo[2] + '\n')
            
            # Create a BytesIO object to store the ZIP file
            zip_buffer = BytesIO()

            # Create a ZipFile object
            with zipfile.ZipFile(
                zip_buffer, "a", zipfile.ZIP_DEFLATED, False
            ) as zip_file:
                # Add the FASTA file to the ZIP file with a custom name
                zip_file.write(fasta_file_name, arcname=fasta_file_name.split("/")[-1])

                # Add the BED file to the ZIP file with a custom name
                zip_file.write(fasta_file_name_oligos, arcname=fasta_file_name_oligos.split("/")[-1])

                # Add the GBK file to the ZIP file with a custom name
                zip_file.write(gbk_file, arcname=gbk_file.split("/")[-1])

            # Move the buffer's position to the beginning to ensure all the data is read
            zip_buffer.seek(0)

            # Return the ZIP file as an attachment
            return send_file(
                zip_buffer,
                download_name=out_dict["files_name"] + ".zip",
                as_attachment=True,
            )         

        # Save selected parameters to input windows
        gene_info_form.text_field.default = out_dict["gene_name"]
        gene_info_form.text_field2.default = out_dict["ncbi_id"]

        gene_info_form.text_field3.default = out_dict["guide_seq"]

        gene_info_form.text_field4.default = out_dict["lha"]
        gene_info_form.text_field5.default = out_dict["rha"]
        gene_info_form.text_field6.default = out_dict["flank"]

        gene_info_form.process()

        cts_info_form.text_field8.default = out_dict['CTS_HA']
        cts_info_form.text_field9.default = out_dict['buffer']
        cts_info_form.text_field10.default = out_dict['scrambled_nt']

        cts_info_form.process()

        date_today = str(date.today())
        save_files_form.text_field7.default = out_dict['files_name']
        save_files_form.process()



        forms = {
            "gene_info_form": gene_info_form,
            "save_files_form":save_files_form,
            "cts_info_form":cts_info_form
        }

        return render_template("home.html", out_dict=out_dict, forms=forms)
    
    gene_info_form.text_field4.default = out_dict["lha"]
    gene_info_form.text_field5.default = out_dict["rha"]
    gene_info_form.text_field6.default = out_dict["flank"]

    gene_info_form.process()

    cts_info_form.text_field8.default = out_dict['CTS_HA']
    cts_info_form.text_field9.default = out_dict['buffer']
    cts_info_form.text_field10.default = out_dict['scrambled_nt']

    cts_info_form.process()

    date_today = str(date.today())
    save_files_form.text_field7.default = out_dict['files_name']
    save_files_form.process()

    return render_template("home.html", out_dict=out_dict, forms=forms)


@app.route("/", methods=["GET", "POST"])
def root():
    """
    Load main page of server
    """
    return index(out_dict)

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)
