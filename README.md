# SubSketch

**SubSketch** generates visualizations of biosynthetic gene clusters (BGCs) and sub-clusters using **ClusterClue** results. It produces both **SVG figures** and **HTML reports** to help explore and present detected **subcluster motifs** across BGCs.

Use SubSketch to create graphics for:
* BGC layouts with CDS features and domain architecture
* Subcluster motif hits with gene-specific highlighting
* Weight matrices of subcluster motifs (positive/negative contributions)
* Chemical structures with substructure matching (RDKit)


## 🚀 Quickstart

```python
import subsketch

# Load everything needed
session = subsketch.SubSketchSession(
    motifs_file="clusterclue_motifs.txt",
    genbank_dir="bgc_gbks",
    domain_hits_file="domain_hits.txt",
    motif_hits_file="detected_motifs.tsv",
    compounds_file="compound_smiles.tsv",
)
session.load()

# Generate HTML report for a BGC
html = session.html_report_for_bgc("BGC0001234")
with open("bgc_report.html", "w") as f:
    f.write(html)

# Generate HTML report for a motif
html = session.html_report_for_motif("1234")
with open("motif_report.html", "w") as f:
    f.write(html)
```

See `examples/` folder for detailed workflows.

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
