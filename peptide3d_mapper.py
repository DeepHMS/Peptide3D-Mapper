import streamlit as st
import pandas as pd
import numpy as np
from Bio import SeqIO
import py3Dmol
import io
import requests
from matplotlib import colormaps
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
import zipfile
from matplotlib import patches
import matplotlib.colors as mcolors  # For rgb2hex

# Set wide layout
st.set_page_config(layout="wide", page_title="Peptide3D Mapper")

def z_score(intensities):
    log_int = np.log10(intensities + 1)
    mean_log = np.mean(log_int)
    std_log = np.std(log_int)
    if std_log == 0:
        return np.zeros_like(log_int)
    else:
        return (log_int - mean_log) / std_log

def map_peptides_to_residues(df, protein_seq, intensity_col, overlap_strategy='merge'):
    if 'Stripped.Sequence' not in df.columns or 'Protein.Group' not in df.columns:
        raise ValueError("CSV must have 'Stripped.Sequence' and 'Protein.Group' columns.")
    
    seq_len = len(protein_seq)
    residue_vals = [None] * seq_len

    peptides = df.groupby('Stripped.Sequence')[intensity_col].mean().reset_index()
    z_scores = z_score(peptides[intensity_col])

    for idx, row in peptides.iterrows():
        pep = row['Stripped.Sequence']
        start = protein_seq.find(pep)
        if start == -1:
            continue
        end = start + len(pep)
        for i in range(start, end):
            z_val = z_scores[idx]
            if residue_vals[i] is None:
                residue_vals[i] = [z_val]
            else:
                residue_vals[i].append(z_val)

    for i in range(seq_len):
        if residue_vals[i]:
            if overlap_strategy == 'merge':
                residue_vals[i] = np.mean(residue_vals[i])
            elif overlap_strategy == 'highest':
                residue_vals[i] = np.max(residue_vals[i])
            elif overlap_strategy in ['none', 'last']:
                residue_vals[i] = residue_vals[i][-1]
            else:
                residue_vals[i] = np.mean(residue_vals[i])  # Default to merge
        else:
            residue_vals[i] = None
    return residue_vals

def generate_colormap(residue_vals, cmap_name='autumn'):
    cmap = colormaps[cmap_name]
    vals = [v for v in residue_vals if v is not None]
    vmin, vmax = (min(vals), max(vals)) if vals else (0, 1)
    hex_colors = []
    for val in residue_vals:
        if val is None:
            hex_colors.append('#d3d3d3')
        else:
            norm = (val - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            rgb = cmap(norm)[:3]
            hex_colors.append(mcolors.rgb2hex(rgb))
    return hex_colors, vmin, vmax

# ✅ Responsive 3D Viewer
def render_viewer(pdb_str, residue_vals, bg_color, title, width_pct="100%", height="420px"):
    hex_colors, vmin, vmax = generate_colormap(residue_vals)
    view = py3Dmol.view()
    view.addModel(pdb_str, 'pdb')
    view.setBackgroundColor(bg_color)
    view.setStyle({}, {'cartoon': {'color': 'lightgray'}})
    for i, c in enumerate(hex_colors):
        view.setStyle({'resi': str(i+1)}, {'cartoon': {'color': c}})
    view.zoomTo()

    st.markdown(f"#### {title}")
    # HTML wrapper for responsiveness
    html_str = f"""
    <div style="width:{width_pct}; height:{height};">
        {view._make_html()}
    </div>
    """
    st.components.v1.html(html_str, height=int(height.replace("px","")))

# ✅ Responsive Linear Plot
def render_linear_plot(residue_vals, title, seq_len, vmin, vmax):
    fig, ax = plt.subplots(figsize=(12, 1))  # Base size, Streamlit scales
    cmap = colormaps['autumn']
    ax.add_patch(patches.Rectangle((0, 0), seq_len, 1, facecolor='lightgray', edgecolor='none'))
    for i in range(seq_len):
        if residue_vals[i] is not None:
            norm = (residue_vals[i] - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            ax.add_patch(patches.Rectangle((i, 0), 1, 1, facecolor=cmap(norm)[:3], edgecolor='none'))
    ax.set_xlim(0, seq_len)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel(f'Amino Acid Position ({title})')
    ax.set_xticks(range(0, seq_len + 1, max(1, seq_len // 10)))
    plt.tight_layout()
    st.pyplot(fig, use_container_width=True)  # ✅ responsive

def create_download_zip(protein_of_interest, pdb_str, peptide_data, residue_data, conditions, min_max_logs, seq_len):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
        # PDB
        zipf.writestr(f"{protein_of_interest}_protein.pdb", pdb_str)

        # Peptide CSVs
        for condition in conditions:
            peptide_csv = peptide_data[condition].to_csv(index=False)
            zipf.writestr(f"{protein_of_interest}_{condition}_peptides.csv", peptide_csv)

        # PyMOL scripts
        cmap = colormaps['autumn']
        for condition in conditions:
            pml_content = f"load {protein_of_interest}_protein.pdb\nhide everything\nshow cartoon\ncolor gray90, all\nzoom\n"
            min_log, max_log = min_max_logs[condition]
            for i in range(seq_len):
                if residue_data[condition][i] is not None:
                    norm = (residue_data[condition][i] - min_log) / (max_log - min_log) if max_log > min_log else 0.5
                    color_hex = mcolors.rgb2hex(cmap(norm)[:3])
                    pml_content += f"color {color_hex}, resi {i+1}\n"
            zipf.writestr(f"{protein_of_interest}_{condition}_pymol_script.pml", pml_content)

        # Linear JPEGs (keep fixed high resolution for downloads)
        for condition in conditions:
            fig, ax = plt.subplots(figsize=(12, 1), dpi=600)
            ax.add_patch(patches.Rectangle((0, 0), seq_len, 1, facecolor='lightgray', edgecolor='none'))
            min_log, max_log = min_max_logs[condition]
            for i in range(seq_len):
                if residue_data[condition][i] is not None:
                    norm = (residue_data[condition][i] - min_log) / (max_log - min_log) if max_log > min_log else 0.5
                    ax.add_patch(patches.Rectangle((i, 0), 1, 1, facecolor=cmap(norm)[:3], edgecolor='none'))
            ax.set_xlim(0, seq_len)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.set_xlabel(f'Amino Acid Position ({condition})')
            ax.set_xticks(range(0, seq_len + 1, max(1, seq_len // 10)))
            img_buffer = io.BytesIO()
            plt.savefig(img_buffer, format='jpeg', dpi=600, bbox_inches='tight')
            img_buffer.seek(0)
            zipf.writestr(f"{protein_of_interest}_{condition}_linear.jpeg", img_buffer.read())
            plt.close(fig)

    zip_buffer.seek(0)
    return zip_buffer
