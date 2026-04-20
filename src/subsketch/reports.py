from subsketch.draw.bgc import draw_bgc, draw_subcluster_hit, draw_annotated_subcluster
from subsketch.draw.motif import plot_subcluster_motif
from subsketch.draw.molecule import draw_compounds, draw_compounds_with_substruct_flexible
from subsketch.loaders import load_bgc
from pathlib import Path
import textwrap  # ADD THIS IMPORT


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
        text: Main title text (can include HTML)
        subtitle: Subtitle text

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
    motif_reports_base_url: str = None,
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
        # Create title with optional link to motif report
        motif_id = motif_hit['motif_id']
        if motif_reports_base_url:
            motif_title = (
                f"Hit for Subcluster Motif "
                f"<a href='{motif_reports_base_url}/{motif_id}.html' "
                f"style='color: #047857; text-decoration: none;' "
                f"onmouseover=\"this.style.textDecoration='underline'\" "
                f"onmouseout=\"this.style.textDecoration='none'\">"
                f"{motif_id}</a>"
            )
        else:
            motif_title = f"Hit for Subcluster Motif {motif_id}"
        
        html_content += _subcluster_title(
            motif_title,
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
    include_motif_plots: bool = True,
    bgc_reports_base_url: str = None,
):
    html_content = _html_head()
    
    if include_title:
        html_content += f"<h1>Subcluster Motif #{motif_id}</h1>\n"

    for motif_hit in motif_hits:
        bgc_id = motif_hit['bgc_id']

        # Create title with optional link to BGC report
        if bgc_reports_base_url:
            bgc_title = (
                f"Subcluster Motif Hit in "
                f"<a href='{bgc_reports_base_url}/{bgc_id}.html' "
                f"style='color: #1e40af; text-decoration: none;' "
                f"onmouseover=\"this.style.textDecoration='underline'\" "
                f"onmouseout=\"this.style.textDecoration='none'\">"
                f"{bgc_id}</a>"
            )
        else:
            bgc_title = f"Subcluster Motif Hit in {bgc_id}"
        
        html_content += _subcluster_title(
            bgc_title,
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


def generate_master_index_html(
    bgc_entries: list,
    motif_entries: list,
    bgc_reports_url: str = "bgc_reports",
    motif_reports_url: str = "motif_reports"
) -> str:
    """Generate a master index.html page with both BGCs and motifs with search.
    
    Args:
        bgc_entries: List of dicts with keys: href, bgc_id, num_motifs, motif_ids
        motif_entries: List of dicts with keys: href, motif_id, num_bgcs, bgc_ids
        bgc_reports_url: Relative URL to BGC reports directory
        motif_reports_url: Relative URL to motif reports directory
    
    Returns:
        str: Complete HTML content for the master index page
    """
    
    total_bgcs = len(bgc_entries)
    total_motifs = len(motif_entries)
    total_hits = sum(entry['num_motifs'] for entry in bgc_entries)
    
    html_content = textwrap.dedent(f"""\
        <!DOCTYPE html>
        <html lang="en">
        <head>
          <meta charset='utf-8'>
          <meta name="viewport" content="width=device-width, initial-scale=1.0">
          <title>SubClusterMotif Hits</title>
          <style>
            * {{ margin: 0; padding: 0; box-sizing: border-box; }}
            
            body {{
              font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
              line-height: 1.6; color: #333;
              background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
              min-height: 100vh; padding: 2rem;
            }}
            
            .container {{
              max-width: 1600px; margin: 0 auto;
              background: rgba(255, 255, 255, 0.95);
              backdrop-filter: blur(10px); border-radius: 20px;
              box-shadow: 0 20px 40px rgba(0,0,0,0.1); overflow: hidden;
            }}
            
            header {{
              background: linear-gradient(135deg, #4f46e5 0%, #7c3aed 100%);
              color: white; padding: 3rem 2.5rem; text-align: center;
            }}
            
            h1 {{ 
              font-size: 3rem; font-weight: 700; margin-bottom: 0.5rem; 
              text-shadow: 0 2px 4px rgba(0,0,0,0.2); 
            }}
            
            .subtitle {{ opacity: 0.9; font-size: 1.2rem; margin-bottom: 1.5rem; }}
            
            .stats {{ 
              display: flex; justify-content: center; gap: 2rem; 
              margin-top: 1.5rem; flex-wrap: wrap; 
            }}
            
            .stat {{ 
              text-align: center; background: rgba(255,255,255,0.2);
              padding: 1.25rem 2rem; border-radius: 12px; 
              backdrop-filter: blur(5px); min-width: 140px;
            }}
            
            .stat-number {{ font-size: 2.5rem; font-weight: 800; display: block; }}
            .stat-label {{ font-size: 0.9rem; opacity: 0.9; text-transform: uppercase; 
                          letter-spacing: 0.05em; margin-top: 0.25rem; }}
            
            main {{ padding: 2.5rem; }}
            
            .search-section {{
              margin-bottom: 2rem;
              background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 100%);
              padding: 2rem;
              border-radius: 16px;
              box-shadow: 0 4px 16px rgba(0,0,0,0.05);
            }}
            
            .search-box {{
              display: flex;
              gap: 1rem;
              margin-bottom: 1rem;
            }}
            
            .search-input {{
              flex: 1;
              padding: 1rem 1.5rem;
              font-size: 1rem;
              border: 2px solid #cbd5e1;
              border-radius: 12px;
              outline: none;
              transition: all 0.2s ease;
            }}
            
            .search-input:focus {{
              border-color: #4f46e5;
              box-shadow: 0 0 0 3px rgba(79, 70, 229, 0.1);
            }}
            
            .search-filters {{
              display: flex;
              gap: 1rem;
              flex-wrap: wrap;
            }}
            
            .filter-btn {{
              padding: 0.5rem 1.25rem;
              border: 2px solid #cbd5e1;
              background: white;
              border-radius: 8px;
              cursor: pointer;
              transition: all 0.2s ease;
              font-weight: 500;
              font-size: 0.9rem;
            }}
            
            .filter-btn.active {{
              background: #4f46e5;
              color: white;
              border-color: #4f46e5;
            }}
            
            .filter-btn:hover:not(.active) {{
              border-color: #4f46e5;
              background: #f8fafc;
            }}
            
            .tabs {{
              display: flex;
              gap: 1rem;
              margin-bottom: 2rem;
              border-bottom: 2px solid #e2e8f0;
            }}
            
            .tab {{
              padding: 1rem 2rem;
              cursor: pointer;
              border: none;
              background: none;
              font-size: 1.1rem;
              font-weight: 600;
              color: #64748b;
              position: relative;
              transition: all 0.2s ease;
            }}
            
            .tab:hover {{
              color: #4f46e5;
            }}
            
            .tab.active {{
              color: #4f46e5;
            }}
            
            .tab.active::after {{
              content: '';
              position: absolute;
              bottom: -2px;
              left: 0;
              right: 0;
              height: 2px;
              background: #4f46e5;
            }}
            
            .tab-content {{
              display: none;
            }}
            
            .tab-content.active {{
              display: block;
            }}
            
            .table-container {{ 
              overflow-x: auto; border-radius: 16px;
              box-shadow: 0 8px 32px rgba(0,0,0,0.1); 
              margin-top: 1rem;
            }}
            
            table {{ 
              width: 100%; border-collapse: collapse; background: white;
              border-radius: 16px; overflow: hidden; 
            }}
            
            th {{ 
              background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 100%);
              padding: 1.25rem 1.5rem; text-align: left; font-weight: 600;
              font-size: 0.95rem; text-transform: uppercase; 
              letter-spacing: 0.05em; color: #475569; 
              border-bottom: 3px solid #e2e8f0;
            }}
            
            td {{ 
              padding: 1.25rem 1.5rem; border-bottom: 1px solid #f1f5f9;
              vertical-align: middle; 
            }}
            
            tr:hover {{ 
              background: linear-gradient(90deg, #f8fafc 0%, #f1f5f9 100%);
            }}
            
            .bgc-link {{ 
              color: #1e40af; text-decoration: none; font-weight: 600;
              padding: 0.5rem 0.75rem; border-radius: 8px;
              transition: all 0.2s ease; display: inline-flex;
              align-items: center; gap: 0.5rem; 
            }}
            
            .bgc-link:hover {{ 
              background: #1e40af; color: white;
              transform: translateY(-1px);
              box-shadow: 0 4px 12px rgba(30, 64, 175, 0.4); 
            }}
            
            .motif-link {{ 
              color: #047857; text-decoration: none; font-weight: 600;
              padding: 0.5rem 0.75rem; border-radius: 8px;
              transition: all 0.2s ease; display: inline-flex;
              align-items: center; gap: 0.5rem; 
            }}
            
            .motif-link:hover {{ 
              background: #047857; color: white;
              transform: translateY(-1px);
              box-shadow: 0 4px 12px rgba(4, 120, 87, 0.4); 
            }}

            .id-list {{ 
              background: #f8fafc; padding: 0.75rem; border-radius: 8px;
              font-family: 'Monaco', 'Menlo', monospace; font-size: 0.85rem;
              border-left: 4px solid #4f46e5; word-break: break-all;
              max-height: 100px; overflow-y: auto;
            }}
            
            .id-list a {{
              color: #4f46e5;
              text-decoration: none;
              margin-right: 0.5rem;
              font-weight: 500;
            }}
            
            .id-list a:hover {{
              text-decoration: underline;
            }}
            
            .badge {{ 
              display: inline-block; background: #4f46e5; color: white;
              padding: 0.25rem 0.75rem; border-radius: 12px; font-size: 0.875rem;
              font-weight: 600; 
            }}
            
            .badge.bgc {{ background: #1e40af; }}
            .badge.motif {{ background: #047857; }}
            
            .no-results {{
              text-align: center;
              padding: 3rem;
              color: #64748b;
              font-size: 1.1rem;
            }}
            
            .no-results::before {{
              content: '🔍';
              display: block;
              font-size: 3rem;
              margin-bottom: 1rem;
            }}
            
            @media (max-width: 768px) {{
              body {{ padding: 1rem; }}
              h1 {{ font-size: 2rem; }}
              main {{ padding: 1.5rem; }}
              th, td {{ padding: 1rem 0.75rem; }}
              .stats {{ gap: 1rem; }}
              .stat {{ min-width: 100px; padding: 1rem 1.5rem; }}
            }}
          </style>
        </head>
        <body>
          <div class="container">
            <header>
              <h1>SubClusterMotif Hits</h1>
              <div class="stats">
                <div class="stat">
                  <span class="stat-number">{total_bgcs}</span>
                  <span class="stat-label">BGCs</span>
                </div>
                <div class="stat">
                  <span class="stat-number">{total_motifs}</span>
                  <span class="stat-label">Motifs</span>
                </div>
                <div class="stat">
                  <span class="stat-number">{total_hits}</span>
                  <span class="stat-label">Hits</span>
                </div>
              </div>
            </header>
            <main>
              <div class="search-section">
                <div class="search-box">
                  <input 
                    type="text" 
                    id="searchInput" 
                    class="search-input" 
                    placeholder="Search by BGC ID, Motif ID, or associated IDs..."
                  >
                </div>
                <div class="search-filters">
                  <button class="filter-btn active" data-filter="all">All</button>
                  <button class="filter-btn" data-filter="bgc">BGCs Only</button>
                  <button class="filter-btn" data-filter="motif">Motifs Only</button>
                </div>
              </div>
              
              <div class="tabs">
                <button class="tab active" data-tab="bgc">BGC Reports</button>
                <button class="tab" data-tab="motif">Motif Reports</button>
              </div>
              
              <div id="bgc-tab" class="tab-content active">
                <div class="table-container">
                  <table id="bgcTable">
                    <thead>
                      <tr>
                        <th>BGC ID</th>
                        <th>Detected Motifs</th>
                        <th>Motif IDs</th>
                      </tr>
                    </thead>
                    <tbody>""")
    
    # BGC table rows
    for entry in bgc_entries:
        motif_links = []
        for mid in entry["motif_ids"]:
            motif_links.append(
                f'<a href="{motif_reports_url}/{mid}.html">{mid}</a>'
            )
        motif_ids_str = " ".join(motif_links) if motif_links else "None"
        
        html_content += textwrap.dedent(f"""\
                      <tr>
                        <td><a href="{bgc_reports_url}/{entry['href']}" class="bgc-link">{entry['bgc_id']}</a></td>
                        <td><span class="badge bgc">{entry['num_motifs']}</span></td>
                        <td><div class="id-list">{motif_ids_str}</div></td>
                      </tr>""")
    
    html_content += textwrap.dedent("""\
                    </tbody>
                  </table>
                </div>
                <div class="no-results" id="bgcNoResults" style="display: none;">
                  No BGCs found matching your search
                </div>
              </div>
              
              <div id="motif-tab" class="tab-content">
                <div class="table-container">
                  <table id="motifTable">
                    <thead>
                      <tr>
                        <th>Motif ID</th>
                        <th>BGCs Found</th>
                        <th>BGC IDs</th>
                      </tr>
                    </thead>
                    <tbody>""")
    
    # Motif table rows
    for entry in motif_entries:
        bgc_links = []
        for bid in entry["bgc_ids"]:
            bgc_links.append(
                f'<a href="{bgc_reports_url}/{bid}.html">{bid}</a>'
            )
        bgc_ids_str = " ".join(bgc_links) if bgc_links else "None"
        
        html_content += textwrap.dedent(f"""\
                      <tr>
                        <td><a href="{motif_reports_url}/{entry['href']}" class="motif-link">{entry['motif_id']}</a></td>
                        <td><span class="badge motif">{entry['num_bgcs']}</span></td>
                        <td><div class="id-list">{bgc_ids_str}</div></td>
                      </tr>""")
    
    html_content += textwrap.dedent("""\
                    </tbody>
                  </table>
                </div>
                <div class="no-results" id="motifNoResults" style="display: none;">
                  No motifs found matching your search
                </div>
              </div>
            </main>
          </div>
          
          <script>
            // Tab switching
            const tabs = document.querySelectorAll('.tab');
            const tabContents = document.querySelectorAll('.tab-content');
            
            tabs.forEach(tab => {
              tab.addEventListener('click', () => {
                const tabName = tab.dataset.tab;
                
                tabs.forEach(t => t.classList.remove('active'));
                tabContents.forEach(tc => tc.classList.remove('active'));
                
                tab.classList.add('active');
                document.getElementById(`${tabName}-tab`).classList.add('active');
              });
            });
            
            // Filter buttons
            const filterBtns = document.querySelectorAll('.filter-btn');
            let currentFilter = 'all';
            
            filterBtns.forEach(btn => {
              btn.addEventListener('click', () => {
                filterBtns.forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
                currentFilter = btn.dataset.filter;
                performSearch();
              });
            });
            
            // Search functionality
            const searchInput = document.getElementById('searchInput');
            const bgcTable = document.getElementById('bgcTable');
            const motifTable = document.getElementById('motifTable');
            const bgcNoResults = document.getElementById('bgcNoResults');
            const motifNoResults = document.getElementById('motifNoResults');
            
            searchInput.addEventListener('input', performSearch);
            
            function performSearch() {{
              const searchTerm = searchInput.value.toLowerCase();
              
              if (currentFilter === 'all' || currentFilter === 'bgc') {{
                filterTable(bgcTable, searchTerm, bgcNoResults);
              }}
              
              if (currentFilter === 'all' || currentFilter === 'motif') {{
                filterTable(motifTable, searchTerm, motifNoResults);
              }}
              
              // Hide/show tables based on filter
              if (currentFilter === 'bgc') {{
                document.getElementById('bgc-tab').style.display = 'block';
                document.getElementById('motif-tab').style.display = 'none';
                tabs[0].click();
              }} else if (currentFilter === 'motif') {{
                document.getElementById('bgc-tab').style.display = 'none';
                document.getElementById('motif-tab').style.display = 'block';
                tabs[1].click();
              }} else {{
                document.getElementById('bgc-tab').style.display = 'block';
                document.getElementById('motif-tab').style.display = 'block';
              }}
            }}
            
            function filterTable(table, searchTerm, noResultsDiv) {{
              const tbody = table.querySelector('tbody');
              const rows = tbody.querySelectorAll('tr');
              let visibleCount = 0;
              
              rows.forEach(row => {{
                const text = row.textContent.toLowerCase();
                if (text.includes(searchTerm)) {{
                  row.style.display = '';
                  visibleCount++;
                }} else {{
                  row.style.display = 'none';
                }}
              }});
              
              // Show/hide no results message
              if (visibleCount === 0) {{
                table.parentElement.style.display = 'none';
                noResultsDiv.style.display = 'block';
              }} else {{
                table.parentElement.style.display = 'block';
                noResultsDiv.style.display = 'none';
              }}
            }}
          </script>
        </body>
        </html>""")
    
    return html_content