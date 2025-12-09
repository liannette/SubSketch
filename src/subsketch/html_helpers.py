from subsketch.bgc import draw_bgc, draw_subcluster_hit, draw_annotated_subcluster
from subsketch.motif import plot_subcluster_motif
from subsketch.loaders import load_bgc, load_mibig_bgc
from subsketch.molecule import draw_compounds, draw_compounds_with_substruct_flexible


def _html_head():
    # change the typeface for the html content
    html_head = (
        "<style>\n"
        "\tbody {\n"
        "\t\tfont-family: 'Arial', sans-serif;\n"
        "\t\tline-height: 1.6;\n"
        "\t}\n"
        "</style>\n"
    )
    return html_head

def _subcluster_title(motif_hit):
    """Wrap SVG content with a title in HTML format.

    Args:
        motif_hit (dict): The motif hit information containing details for the title.

    Returns:
        str: HTML string containing the title and SVG content.
    """
    title = f"{motif_hit['bgc_id']} - Subcluster Motif #{motif_hit['motif_id']}"
    subtitle = f"Threshold: {motif_hit['threshold']}\nScore: {motif_hit['score']}"

    html_content = (
        "<div style=\"margin: 1.5em 0 1em 0;\">\n"
        "    <h3 style=\"margin: 0 0 0.3em 0; line-height: 1.2;\">\n"
        f"        {title}\n"
        "        <br><small style=\"font-size: 0.8em; color: #666; font-weight: 400;\">\n"
        f"            {subtitle}\n"
        "        </small>\n"
        "    </h3>\n"
        "</div>\n"
    )
    return html_content


def generate_html_report_for_bgc(
    bgc_id: str,
    bgc_length: int,
    cds_features: list,
    bgc_detected_motifs: list,
    bgc_domains_hits: dict,
    domain_colors: dict,
    compounds: list = [],
    motifs: dict = {},
    include_compound_plots: bool = True,
    include_motif_plots: bool = True,
    include_bgc_plot: bool = True,
    scaling: float = 30,
):
    html_content = _html_head()

    # Draw BGC
    html_content += f"<h1>{bgc_id}</h1>\n"
    if include_bgc_plot:
        bgc_svg = draw_bgc(
            bgc_length=bgc_length,
            cds_features=cds_features,
            color_genes=True,
            color_domains=False,
            scaling=scaling,
        )
        html_content += f"<div>{bgc_svg}</div>"
    
    if include_compound_plots:
        # Draw BGC compounds
        compounds_svg = draw_compounds(compounds)
        html_content += f"<div>{compounds_svg}</div>"

    # Draw subcluster hits
    for motif_hit in bgc_detected_motifs:
        html_content += _subcluster_title(motif_hit)

        subcluster_svg = draw_subcluster_hit(
            bgc_length=bgc_length,
            cds_features=cds_features,
            motif_hit=motif_hit,
            bgc_domain_hits=bgc_domains_hits,
            domain_colors=domain_colors,
            scaling=scaling,
        )
        html_content += f"<div class='subcluster_hit'>{subcluster_svg}</div>"

        if include_motif_plots:
            motif_svg = plot_subcluster_motif(
                motifs[motif_hit['motif_id']],
                motif_hit,
            )
            html_content += f"<div class='motif_plot'>{motif_svg}</div>"

    return html_content


def generate_html_for_motif(
    motif_id,
    motif_hits, 
    gbks_dirpath, 
    domains, 
    domain_colors, 
    motifs, 
    compounds_info=None, 
    scaling=30, 
    include_motif_plots=True
):
    
    html_content = _html_head()
    html_content += f"<h1>Subcluster Motif #{motif_id}</h1>\n"

    for motif_hit in motif_hits:

        bgc_id = motif_hit['bgc_id']
        html_content += _subcluster_title(motif_hit)

        if compounds_info:
            bgc_compounds = compounds_info.get(bgc_id, [])
            compounds_svg = draw_compounds(bgc_compounds)
            html_content += f"<div>{compounds_svg}</div>"

        bgc_data = load_bgc(gbks_dirpath / f"{bgc_id}.gbk")
        bgc_domains = domains.get(bgc_id, {})   
        subcluster_svg = draw_subcluster_hit(
            bgc_data=bgc_data,
            motif_hit=motif_hit,
            bgc_domain_hits=bgc_domains,
            domain_colors=domain_colors,
            scaling=scaling,
        )
        html_content += f"<div class='subcluster_hit'>{subcluster_svg}</div>"
        
        if include_motif_plots:
            motif_svg = plot_subcluster_motif(
                motifs[motif_hit['motif_id']],
                motif_hit,
            )

        html_content += f"<div class='motif_plot'>{motif_svg}</div>"

    return html_content


def generate_html_for_subcluster(subcluster, gbks_dirpath, compounds, scaling=30):

    subcluster_id = subcluster["subcluster_id"]
    bgc_id = subcluster["bgc_id"]
    substructure_name = subcluster["substructure_name"]
    substructure_smiles = subcluster["substructure_smiles"]
    original_sequence = subcluster.get("orig_seq", None)

    bgc_data = load_mibig_bgc(gbks_dirpath / f"{bgc_id}.gbk")

    # Generate HTML content
    html_content = _html_head()

    # Title and subtitle
    title = f"SC{subcluster_id.zfill(6)}"
    subtitle = f"{bgc_id}: {bgc_data['description']} from {bgc_data['organism'] }"
    link = f"https://mibig.secondarymetabolites.org/repository/{bgc_id}"
    html_content += (
        "<div style=\"margin: 1.5em 0 1em 0;\">\n"
        "    <h1 style=\"margin: 0 0 0.3em 0; line-height: 1.2;\">\n"
        f"        {title}\n"
        "        <br>\n"
        "        <small style=\"font-size: 0.5em; color: #666; font-weight: 400;\">\n"
        "        From: \n"
        f"            <a href=\"{link}\" style=\"color: #007bff; text-decoration: none;\"\n"
        f"            onmouseover=\"this.style.textDecoration='underline'\"\n"
        f"            onmouseout=\"this.style.textDecoration='none'\">\n"
        f"                {subtitle}\n"
        "            </a>\n"
        "        </small>\n"
        "    </h1>\n"
        "</div>\n"
    )

    # Draw BGC compounds with highlighted substructure (only first compound)
    compounds_highlighted_svg = draw_compounds_with_substruct_flexible(
        compounds=compounds[:1],
        substruct_smiles=subcluster["substructure_smiles"],
        show_names=False,
        scaling=2,
    )
    html_content += f"<div style='text-align: center;'>{compounds_highlighted_svg}</div>"

   # Draw subcluster in BGC context
    subcluster_svg = draw_annotated_subcluster(
        annotated_subcluster=subcluster,
        bgc_data=bgc_data,
        scaling=scaling,
    )
    html_content += f"<div style='text-align: center;'>{subcluster_svg}</div>"

    # Create a table for information
    genbank_accession = original_sequence
    genes_with_protein_ids = [f"{pid} ({gene})" for gene, pid in zip(subcluster['genes'], subcluster['protein_ids'])]
    genes_html = "<br>".join(str(gene) for gene in genes_with_protein_ids)
    substructure_name = subcluster["substructure_name"]
    substructure_svg = draw_compounds([(substructure_name, substructure_smiles)], show_names=False, scaling=0.8)
    substructure_class = subcluster.get("substructure_class", "N/A")
    pathway_quality = subcluster.get("pathway_quality", "N/A")
    references = ["PMID: " + pmid for pmid in subcluster.get("pubmed_id", [])]
    references_html = '<br>'.join(references)


    html_content += (
        "<table style='width: 80%; margin: 1em auto; border-collapse: collapse;'>"
        f"<tr>"
        "<td style='border: 1px solid #ddd; padding: 8px;'>Protein IDs (Genes)</td>"
        f"<td style='border: 1px solid #ddd; padding: 8px; font-family: monospace; white-space: pre-wrap;'>{genes_html}</td>"
        "</tr>"
        f"<tr>"
        "<td style='border: 1px solid #ddd; padding: 8px;'>NCBI GenBank Sequence Accession</td>"
        f"<td style='border: 1px solid #ddd; padding: 8px; font-family: monospace;'>{genbank_accession}</td>"
        "</tr>"
        f"<tr>"
        "<td style='border: 1px solid #ddd; padding: 8px;'>Associated Substructure</td>"
        f"<td style='border: 1px solid #ddd; padding: 8px;'>{substructure_svg}<br>{substructure_name}</td>"
        "</tr>"
        f"<tr>"
        "<td style='border: 1px solid #ddd; padding: 8px;'>Substructure Class</td>"
        f"<td style='border: 1px solid #ddd; padding: 8px;'>{substructure_class}</td>"
        "</tr>"
        f"<tr>"
        "<td style='border: 1px solid #ddd; padding: 8px;'>Pathway Annotation Quality</td>"
        f"<td style='border: 1px solid #ddd; padding: 8px;'>{pathway_quality}</td>"
        "</tr>"
        f"<tr>"
        "<td style='border: 1px solid #ddd; padding: 8px;'>References</td>"
        f"<td style='border: 1px solid #ddd; padding: 8px; font-family: monospace;'>{references_html}</td>"
        "</tr>"

        "</table>"
    )


    return html_content