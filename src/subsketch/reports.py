from subsketch.draw.bgc import draw_bgc, draw_subcluster_hit, draw_annotated_subcluster
from subsketch.draw.motif import plot_subcluster_motif
from subsketch.draw.molecule import draw_compounds, draw_compounds_with_substruct_flexible
from subsketch.loaders import load_bgc
from pathlib import Path


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


def _subcluster_title(text: str, subtitle: str):
    """Wrap SVG content with a title in HTML format.

    Args:
        motif_hit (dict): The motif hit information containing details for the title.

    Returns:
        str: HTML string containing the title and SVG content.
    """
    html_content = (
        "<div style=\"margin: 1.5em 0 1em 0;\">\n"
        "    <h3 style=\"margin: 0 0 0.3em 0; line-height: 1.2;\">\n"
        f"        {text}\n"
        f"        <br><small style=\"font-size: 0.8em; color: #666; font-weight: 400;\">\n"
        f"            {subtitle}\n"
        "        </small>\n"
        "    </h3>\n"
        "</div>\n"
    )
    return html_content


# def _hit_info(text):
#     """Wrap SVG content with a title in HTML format.

#     Args:
#         motif_hit (dict): The motif hit information containing details for the title.

#     Returns:
#         str: HTML string containing the title and SVG content.
#     """
#     html_content = (
#         "<div style=\"margin: 1.5em 0 1em 0;\">\n"
#         "    <small style=\"font-size: 0.8em; color: #666; font-weight: 400;\">\n"
#         f"        {text}\n"
#         "    </small>\n"
#         "</div>\n"
#     )
#     return html_content

# def _subcluster_title_for_motif_report(motif_hit):
#     """Wrap SVG content with a title in HTML format.

#     Args:
#         motif_hit (dict): The motif hit information containing details for the title.

#     Returns:
#         str: HTML string containing the title and SVG content.
#     """
#     title = f"Hit in {motif_hit['bgc_id']}"
#     subtitle = f"Threshold: {motif_hit['threshold']}\nScore: {motif_hit['score']}"

#     html_content = (
#         "<div style=\"margin: 1.5em 0 1em 0;\">\n"
#         "    <h3 style=\"margin: 0 0 0.3em 0; line-height: 1.2;\">\n"
#         f"        {title}\n"
#         "        <br><small style=\"font-size: 0.8em; color: #666; font-weight: 400;\">\n"
#         f"            {subtitle}\n"
#         "        </small>\n"
#         "    </h3>\n"
#         "</div>\n"
#     )
#     return html_content


def toggleable_panel_html(panel_content: str, main_content: str, panel_id: str) -> str:
    html_content = (
        "<style>\n"
        "  .toggle-panel {\n"
        "    display: none;\n"
        "  }\n"
        "  .panel {\n"
        "    max-height: 0;\n"
        "    overflow: hidden;\n"
        "    transition: max-height 0.5s ease-out;\n"
        "    background-color: #f1f1f1;\n"
        "    padding: 0 18px;\n"
        "  }\n"
        "  .toggle-panel:checked + label + .panel {\n"
        "    max-height: 500px;\n"
        "    padding: 18px;\n"
        "  }\n"
        "  .toggle-btn {\n"
        "    cursor: pointer;\n"
        "    color: #003366;\n"
        "    text-decoration: underline;\n"
        "    font-size: 0.8em;\n"
        "  }\n"
        "</style>\n"
        f"<input type=\"checkbox\" class=\"toggle-panel\" id=\"{panel_id}\">\n"
        f"<label for=\"{panel_id}\" class=\"toggle-btn\">{main_content}</label>\n"
        f"<div class=\"panel\">\n{panel_content}\n</div>\n"
    )
    return html_content


def generate_html_report_for_bgc(
    bgc_data: dict,
    bgc_detected_motifs: list,
    bgc_domains_hits: dict,
    domain_colors: dict,
    compounds: list = [],
    motifs: dict = {},
    include_title: bool = True,
    include_compound_plots: bool = True,
    include_motif_plots: bool = True,
    include_bgc_plot: bool = True,
    scaling: float = 30,
):
    html_content = _html_head()

    if include_title:
        html_content += f"<h1>{bgc_data['id']}</h1>\n"

    # Draw BGC compounds
    if compounds and include_compound_plots:
        compounds_svg = draw_compounds(compounds, scaling=1.5)
        html_content += f"<div>{compounds_svg}</div>"

    # Draw BGC
    if include_bgc_plot:
        bgc_svg = draw_bgc(
            bgc_data=bgc_data,
            color_genes=True,
            color_domains=False,
            scaling=scaling,
        )
        html_content += f"<div>{bgc_svg}</div>"

    # Draw subcluster hits
    for motif_hit in bgc_detected_motifs:
        html_content += _subcluster_title(
            f"Hit for Subcluster Motif {motif_hit['motif_id']}", 
            f"Score: {motif_hit['score']} | Threshold: {motif_hit['threshold']}"
            )

        subcluster_svg = draw_subcluster_hit(
            motif_hit=motif_hit,
            bgc_data=bgc_data,
            bgc_domain_hits=bgc_domains_hits,
            domain_colors=domain_colors,
            scaling=scaling,
        )
        html_content += f"<div class='subcluster_hit'>{subcluster_svg}</div>"
        # html_content += _hit_info(f"Score: {motif_hit['score']} | Threshold: {motif_hit['threshold']}")

        if motifs and include_motif_plots:
            motif_svg = plot_subcluster_motif(
                motifs[motif_hit['motif_id']],
                motif_hit,
            )
            html_content += toggleable_panel_html(
                panel_content=motif_svg,
                main_content="<b>View Motif Plot</b>",
                panel_id=f"panel_{bgc_data['id']}_{motif_hit['motif_id']}"
            )
            #f"<div class='motif_plot'>{motif_svg}</div>"

    return html_content


def generate_html_for_motif(
    motif_id: str,
    motif_hits: list,
    gbks_dirpath: str | Path,
    domains: dict,
    domain_colors: dict,
    motifs: dict = {},
    compounds_info: dict = {}, 
    scaling: int = 30, 
    include_title: bool = True,
    include_compound_plots: bool = True,
    include_motif_plots: bool = True
):
    html_content = _html_head()
    
    if include_title:
        html_content += f"<h1>Subcluster Motif #{motif_id}</h1>\n"

    for motif_hit in motif_hits:

        bgc_id = motif_hit['bgc_id']
        html_content += _subcluster_title(
            f"Subcluster Motif Hit in {bgc_id}",
            f"Score: {motif_hit['score']} | Threshold: {motif_hit['threshold']}"
            )

        if compounds_info and include_compound_plots:
            bgc_compounds = compounds_info.get(bgc_id, [])
            if bgc_compounds:
                compounds_svg = draw_compounds(bgc_compounds)
                html_content += f"<div>{compounds_svg}</div>"

        bgc_data = load_bgc(Path(gbks_dirpath) / f"{bgc_id}.gbk")
        bgc_domains = domains.get(bgc_id, {})   
        subcluster_svg = draw_subcluster_hit(
            bgc_data=bgc_data,
            motif_hit=motif_hit,
            bgc_domain_hits=bgc_domains,
            domain_colors=domain_colors,
            scaling=scaling,
        )
        html_content += f"<div class='subcluster_hit'>{subcluster_svg}</div>"
        
        if motifs and include_motif_plots:
            motif_svg = plot_subcluster_motif(
                motifs[motif_hit['motif_id']],
                motif_hit,
            )

            html_content += toggleable_panel_html(
                panel_content=motif_svg,
                main_content="<b>View Motif Plot</b>",
                panel_id=f"panel_{motif_hit['motif_id']}_{bgc_id}"
            )
    return html_content


def generate_html_for_annotated_subcluster(subcluster, bgc_data, compounds, gene_arrow_scaling=30):

    subcluster_id = str(subcluster["id"])
    bgc_id = subcluster["bgc_id"]
    substructure_name = subcluster["substructure_name"]
    substructure_smiles = subcluster["substructure_smiles"]
    original_sequence = subcluster.get("orig_seq", None)

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
        scaling=gene_arrow_scaling,
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