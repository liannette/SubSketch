# SubSketch

**SubSketch** generates high-quality, publication-ready visualizations of biosynthetic gene clusters (BGCs) and sub-clusters using **ClusterClue** results. It produces both **SVG figures** and **HTML reports** to help you explore and present **subcluster motifs** and relationships across BGCs.

Use SubSketch to create clear, customizable graphics that highlight:
* BGC layouts with CDS features and domain architecture
* Subcluster motif hits with gene-specific highlighting
* Weight matrices of subcluster motifs (positive/negative contributions)
* Chemical structures with substructure matching (RDKit)


## 🚀 Quickstart

```python
import subsketch 

# Load everything needed
bgc_data = subsk.loaders.load_bgc("BGC0001234.gbk")
bgc2hits, _ = subsketch.io.read_detected_motifs("detected_motifs.tsv")
domains = subsketch.io.read_domain_hits("domain_hits.txt")
compounds = subsketch.io.read_compounds("compound_smiles.tsv")
motifs = subsketch.io.read_motifs(motifs_filepath)
domain_colors = subsketch.io.read_color_domains_file()

# Generate HTML report
html = generate_html_report_for_bgc(
    bgc_data=bgc_data,
    bgc_domains_hits=domains.get(bgc_data["id"], {}),
    bgc_detected_motifs=bgc2hits.get(bgc_data["id"], []),
    compounds=compounds.get(bgc_data["id"], []),
    motifs=motifs,
    domain_colors=domain_colors,
)

# Save to file
with open("bgc_report.html", "w") as f:
    f.write(html)
```

See `notebook/` folder for detailed workflows.

## 📦 Installation

```bash
# Recommended: Create dedicated environment
conda create -n subsketch python=3.11 -y
conda activate subsketch

# Clone the git repo
git clone https://github.com/liannette/SubSketch.git
cd SubSketch

# Install SubSketch
pip install .
```

To run SubSketch in a jupyter notebook, also install:
```bash
conda install jupyter ipykernel -y
python -m ipykernel install --user --name=subsketch
```
