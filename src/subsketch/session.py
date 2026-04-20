from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Any, Optional
from importlib.resources import files
from subsketch import io, loaders
from subsketch.reports import generate_html_report_for_bgc, generate_html_for_motif


@dataclass
class GlobalData:
    """All data loaded once for a dataset."""
    bgcs: Dict[str, dict]
    domain_hits: Dict[str, Dict[int, List[dict]]]
    bgc2hits: Dict[str, List[dict]]
    motif2hits: Dict[str, List[dict]]
    motifs: Dict[str, Any]
    compounds: Dict[str, List[tuple]]  # optional
    domain_colors: Dict[str, List[int]]


class SubSketchSession:
    """
    High-level interface for SubSketch.

    Load all input files once, then generate HTML reports
    for individual BGCs by ID.
    """

    def __init__(
        self,
        genbank_dir: str | Path,
        domain_hits_file: str | Path,
        motif_hits_file: str | Path,
        motifs_file: str | Path,
        compounds_file: Optional[str | Path] = None,
        domain_colors_file: Optional[str | Path] = None,
    ):
        self.genbank_dir = Path(genbank_dir)
        self.domain_hits_file = Path(domain_hits_file)
        self.motif_hits_file = Path(motif_hits_file)
        self.motifs_file = Path(motifs_file)
        self.compounds_file = Path(compounds_file) if compounds_file else None
        if domain_colors_file:
            self.domain_colors_file = Path(domain_colors_file) 
        else:
            self.domain_colors_file = Path(files("subsketch").joinpath("data")).joinpath("domain_colors.txt")

        self.data: Optional[GlobalData] = None

    def load(self) -> None:
        """Load all data."""
        # BGCs: {bgc_id: {"id", "cds_features", "length", "record", ...}}
        bgcs: Dict[str, dict] = {}
        for gbk_path in sorted(self.genbank_dir.glob("*.gbk")):
            bgc = loaders.load_bgc(gbk_path)
            bgcs[bgc["id"]] = bgc

        domain_hits = io.read_domain_hits(self.domain_hits_file)
        bgc2hits, motif2hits = io.read_detected_motifs(self.motif_hits_file)
        motifs = io.read_motifs(self.motifs_file)
        compounds = io.read_compounds(self.compounds_file) if self.compounds_file else {}

        domain_colors = loaders.load_domain_colors(self.domain_colors_file)
        new_colors = loaders.new_domain_colors(domain_hits, domain_colors)
        if new_colors:
            domain_colors.update(new_colors)
            io.write_domain_colors(domain_colors, self.domain_colors_file)

        self.data = GlobalData(
            bgcs=bgcs,
            domain_hits=domain_hits,
            bgc2hits=bgc2hits,
            motif2hits=motif2hits,
            motifs=motifs,
            compounds=compounds,
            domain_colors=domain_colors,
        )

    def list_bgcs(self) -> List[str]:
        """Return all BGC IDs available in this session."""
        if self.data is None:
            raise RuntimeError("Session not loaded. Call .load() first.")
        return sorted(self.data.bgcs.keys())

    def html_report_for_bgc(
        self,
        bgc_id: str,
        gene_arrow_scaling: int = 30,
        include_title: bool = True,
        include_bgc_plot: bool = True,
        include_compound_plots: bool = True,
        include_motif_plots: bool = True,
        motif_reports_base_url: str = None,
    ) -> str:
        """Generate an HTML report for a single BGC ID."""
        if self.data is None:
            raise RuntimeError("Session not loaded. Call .load() first.")
        if bgc_id not in self.data.bgcs:
            raise KeyError(f"BGC ID not found in session: {bgc_id}")

        bgc_data = self.data.bgcs[bgc_id]
        motif_hits = self.data.bgc2hits.get(bgc_id, [])

        return generate_html_report_for_bgc(
            bgc_data=bgc_data,
            bgc_domains_hits=self.data.domain_hits.get(bgc_id, {}),
            bgc_detected_motifs=motif_hits,
            compounds=self.data.compounds.get(bgc_id, []),
            motifs=self.data.motifs,
            domain_colors=self.data.domain_colors,
            scaling=gene_arrow_scaling,
            include_title=include_title,
            include_motif_plots=include_motif_plots,
            include_bgc_plot=include_bgc_plot,
            include_compound_plots=include_compound_plots,
            motif_reports_base_url=motif_reports_base_url,
        )
    
    def html_report_for_motif(
        self,
        motif_id: str,
        gene_arrow_scaling: int = 30,
        include_title: bool = True,
        include_compound_plots: bool = True,
        include_motif_plots: bool = True,
        bgc_reports_base_url: str = None,
    ) -> str:
        """Generate an HTML report for a single motif ID."""
        if self.data is None:
            raise RuntimeError("Session not loaded. Call .load() first.")
        if motif_id not in self.data.motifs:
            raise KeyError(f"Motif ID not found in session: {motif_id}")

        motif_hits = self.data.motif2hits.get(motif_id, [])

        return generate_html_for_motif(
            motif_id=motif_id,
            motif_hits=motif_hits,
            gbks_dirpath=self.genbank_dir,
            domains=self.data.domain_hits,
            domain_colors=self.data.domain_colors,
            motifs=self.data.motifs,
            compounds_info=self.data.compounds,
            scaling=gene_arrow_scaling,
            include_title=include_title,
            include_motif_plots=include_motif_plots,
            include_compound_plots=include_compound_plots,
            bgc_reports_base_url=bgc_reports_base_url,
        )

    def generate_all_reports_with_master_index(
        self,
        output_dir: str | Path,
        gene_arrow_scaling: int = 60,
        include_compound_plots: bool = True,
        include_motif_plots: bool = True,
    ) -> None:
        """Generate all reports (BGC, Motif) with a master index page.
        
        Args:
            output_dir: Base directory to write all reports
            gene_arrow_scaling: Scaling factor for gene arrows
            include_compound_plots: Whether to include compound visualizations
            include_motif_plots: Whether to include motif plots
        """
        from subsketch.reports import generate_master_index_html
        
        if self.data is None:
            raise RuntimeError("Session not loaded. Call .load() first.")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Define subdirectories
        bgc_dir = output_dir / "bgc_reports"
        motif_dir = output_dir / "motif_reports"
        
        # Generate BGC reports
        bgc_dir.mkdir(parents=True, exist_ok=True)
        bgc_index_entries = []
        
        for idx, bgc_id in enumerate(self.data.bgcs.keys(), 1):
            html_content = self.html_report_for_bgc(
                bgc_id=bgc_id,
                gene_arrow_scaling=gene_arrow_scaling,
                include_compound_plots=include_compound_plots,
                include_motif_plots=include_motif_plots,
                motif_reports_base_url="../motif_reports",
            )
            
            output_file = bgc_dir / f"{bgc_id}.html"
            with open(output_file, "w") as f:
                f.write(html_content)
            
            motif_hits = self.data.bgc2hits.get(bgc_id, [])
            motif_ids = [hit["motif_id"] for hit in motif_hits]
            
            bgc_index_entries.append({
                "href": f"{bgc_id}.html",
                "bgc_id": bgc_id,
                "num_motifs": len(motif_ids),
                "motif_ids": motif_ids
            })
        
        # Generate motif reports
        motif_dir.mkdir(parents=True, exist_ok=True)
        motif_index_entries = []
        
        for idx, motif_id in enumerate(self.data.motifs.keys(), 1):
            html_content = self.html_report_for_motif(
                motif_id=motif_id,
                gene_arrow_scaling=gene_arrow_scaling,
                include_compound_plots=include_compound_plots,
                include_motif_plots=include_motif_plots,
                bgc_reports_base_url="../bgc_reports",
            )
            
            output_file = motif_dir / f"{motif_id}.html"
            with open(output_file, "w") as f:
                f.write(html_content)
            
            motif_hits = self.data.motif2hits.get(motif_id, [])
            bgc_ids = [hit["bgc_id"] for hit in motif_hits]
            
            motif_index_entries.append({
                "href": f"{motif_id}.html",
                "motif_id": motif_id,
                "num_bgcs": len(bgc_ids),
                "bgc_ids": bgc_ids
            })
            
        # Generate master index
        master_index_html = generate_master_index_html(
            bgc_entries=bgc_index_entries,
            motif_entries=motif_index_entries,
            bgc_reports_url="bgc_reports",
            motif_reports_url="motif_reports"
        )
        
        master_index_path = output_dir / "index.html"
        with open(master_index_path, "w") as f:
            f.write(master_index_html)
