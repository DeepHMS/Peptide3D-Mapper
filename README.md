Peptide3D Mapper is a simple, robust web-based tool designed to map peptide data from proteomics experiments to 3D protein structures, enabling visualization of peptide intensity profiles for deeper insights into protein behavior. Built with Jupyter notebooks and deployed using Voilà, it integrates CSV and FASTA inputs with AlphaFold-derived PDB files, offering interactive features like scaling (log/z-score), overlap strategies, and customizable 2D/3D views. Ideal for bioinformaticians and structural biologists, it provides downloadable outputs (PDB, PML, CSV, JPEG) and is hosted on GitHub with Binder for easy access.
Features

Peptide Mapping: Maps peptide sequences from a CSV file to protein structures using a FASTA reference.

3D Visualization: Renders 3D protein structures with color-coded peptide intensities using py3Dmol.

Interactive Widgets: Adjust species, scaling (log/z-score), overlap strategy (none/merge/highest), isoform combination, and background color.

Downloads: Export results as PDB, PyMOL script (PML), CSV, and JPEG files.

Cloud Access: Launch via Binder without local installation.

Requirements

Python 3.9

Dependencies: py3dmol, biopython, requests, ipywidgets, matplotlib, pandas, scipy

Installation & Setup

Local Development

Clone the Repository:
bashgit clone https://github.com/yourusername/Peptide3D-Mapper.git
cd Peptide3D-Mapper

Install Dependencies:
bashpip install -r requirements.txt

Run Locally:
bashvoila peptide_analysis.ipynb --port=8866
Open your browser at http://localhost:8866.

Cloud Deployment (Binder)

Click the Binder link below to launch the app instantly in your browser:

Launch Peptide3D Mapper

No installation needed; sessions are temporary (up to 2 hours).

Usage

Input Files

Peptide_Demo_File.csv: Contains peptide data (e.g., Protein.Group, Stripped.Sequence, 10_Cells_Intensity).

uniprotkb_Human.fasta: Reference protein sequences (e.g., UniProt IDs like O00560).

Steps

Upload Files:

Use the "Upload CSV" and "Upload FASTA" widgets to provide your data files.


Configure Settings:

Set the species (e.g., "Human").

Search or select a protein (e.g., O00560) using the dropdown.

Choose scaling (log/z-score), overlap strategy (none/merge/highest), and background color.

Optionally combine isoforms or specify custom ones.


Visualize:

View the linear intensity plot and 3D structure with color-coded peptides.
A single colorbar shows the scaled intensity range.


Download:

Click "Download Files" to get PDB, PML, CSV, and JPEG outputs.


Notes

If the 3D structure fails to load (e.g., AlphaFold fetch error), a fallback test structure is used.
Large FASTA files (e.g., 29.3 MB) may slow Binder builds; consider using a subset.

Project Structure
textPeptide3D-Mapper/
├── peptide_analysis.ipynb  # Main notebook with Voilà app
├── requirements.txt        # Python dependencies
├── runtime.txt             # Python version
├── launch_app.py           # Script to open Binder link
└── README.md               # This file

Contributing

Fork the repository.
Create a branch for your feature: git checkout -b feature-name.
Commit changes: git commit -m "Add feature-name".
Push and open a pull request.

License
This project is under the MIT License - see the LICENSE file for details (add a LICENSE file if not present).

Support

Issues: Report bugs or suggest features on the GitHub Issues page.

Acknowledgments

Utilizes py3Dmol for 3D rendering and AlphaFold for PDB data.


Version

Developmental version Release: January 24, 2024

Beta Version Release: September 27, 2025




Character Count: Approximately 1,200 characters, providing detailed guidance without overwhelming users.
Customization: Replace yourusername with your actual GitHub username. Add a LICENSE file if you plan to open-source it.
Binder Link: The link is placeholder; generate it via mybinder.org after pushing the repo.
