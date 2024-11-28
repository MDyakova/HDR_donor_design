# HDR Donor Template Design Tool

A **Flask-based web application** packaged in a Docker container for creating HDR (Homology-Directed Repair) donor templates. This tool processes gene and guide RNA information to generate homology arms, donor sequences, and customizable insert sequences for genome editing experiments.

---

## Key Features

- **Gene Information Integration**: Fetches gene data from Ensemble and NCBI databases.
- **Guide RNA Handling**: Processes guide RNA sequences, calculates insertion positions, and integrates related features.
- **Custom Sequence Editing**: Supports adding, removing, or modifying custom elements like promoters and coding sequences.
- **Output Formats**:
  - Generates donor templates in FASTA, BED, and GBK formats.
  - Bundles outputs into a ZIP archive for easy download.
- **Interactive Forms**: User-friendly forms for specifying input parameters and configurations.
- **Dockerized Deployment**: Runs seamlessly in a Docker container, ensuring consistent setup and scalability.

---

## Requirements

- **Docker**: Install Docker from the [official website](https://docs.docker.com/get-docker/).

---

## Getting Started

### 1. Pull the Docker Image

```bash
docker pull mdyakova/hdr_donor_for_crisprcas9:flaskv1```

### 3. Run the Docker Container

```bash
docker run -p 5000:5000 -v ${PWD}\outputs_flask:/src/static/outputs -d dockerfile:vflask1 
```

- **Port Mapping**: Maps the container’s port 5000 to your local port 5000.
- **Volume Mounts**:
  - `src/static/data`: Stores the Excel file with predefined sequences.
  - `src/static/outputs`: Saves generated donor templates and other output files.

### 4. Access the Web Interface

Open your browser and navigate to:
```
http://localhost:5000
```

---

## Usage

1. **Input Gene Information**:
   - Provide Ensemble gene name, NCBI ID, and guide sequence.
   - Specify the sizes for left and right homology arms (LHA/RHA) and flanking sequences.

2. **Generate Donor Template**:
   - Click the **"Generate Sequence"** button to create the donor template.
   - View or modify the generated sequence.

3. **Add Custom Elements**:
   - Upload FASTA files for custom elements like promoters.
   - Specify options for in-frame insertion or reverse sequence orientation.

4. **Export Files**:
   - Download donor sequences in FASTA, BED, or GBK formats.
   - Retrieve oligo sequences and annotations in a ZIP archive.

---

## Environment Variables (Optional)

To customize the application behavior, set the following environment variables in the Docker container:

- `FLASK_APP`: Defaults to `app.py`.
- `FLASK_ENV`: Set to `development` for debugging or `production` for deployment.
- `SECRET_KEY`: Flask secret key for session management.

## Notes

- Ensure the `all_sequences.xlsx` file is formatted correctly with columns for "Elements," "Names," and "Sequences."
- Update the `config.json` file to adjust colors, default parameters, or initial settings.
- The container runs on port `5000` by default. Change the port mapping in the `docker run` command if necessary.

---

Happy designing! 😊

--- 

Let me know if you need additional sections!


## make docker container
docker build -t dockerfile:vflask1 ./  
docker run -p 5000:5000 -v ${PWD}\outputs_flask:/src/static/outputs -d dockerfile:vflask1   

## docker hub
docker tag dockerfile:vflask1 mdyakova/hdr_donor_for_crisprcas9:flaskv1 
docker login
docker push mdyakova/hdr_donor_for_crisprcas9:flaskv1

## Check code quality
python3 -m black main.py
python3 lint.py --threshold 7

## start unit-tests
python3 -m pytest unit_tests.py