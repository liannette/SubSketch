from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Any, Optional
from importlib.resources import files
from multiprocessing import Pool
import logging

from subsketch import io, loaders
from subsketch.reports import generate_html_report_for_bgc, generate_html_for_motif
from subsketch.reports import generate_master_index_html


logger = logging.getLogger(__name__)


@dataclass
class GlobalData:
    """All data loaded once for a dataset."""
    domain_hits: Dict[str, Dict[int, List[dict]]]
    bgc2hits: Dict[str, List[dict]]
    motif2hits: Dict[str, List[dict]]
    motifs: Dict[str, Any]
    compounds: Dict[str, List[tuple]]
    domain_colors: Dict[str, List[int]]

def _generate_single_bgc_html(args):
    """Worker function to generate a single BGC HTML report."""
    (bgc_id, bgc_data, global_data, bgc_dir, gene_arrow_scaling, 
     include_compound_plots, include_motif_plots) = args
    
    motif_hits = global_data.bgc2hits.get(bgc_id, [])
    
    html_content = generate_html_report_for_bgc(
        bgc_data=bgc_data,
        bgc_domains_hits=global_data.domain_hits.get(bgc_id, {}),
        bgc_detected_motifs=motif_hits,
        compounds=global_data.compounds.get(bgc_id, []),
        motifs=global_data.motifs,
        domain_colors=global_data.domain_colors,
        scaling=gene_arrow_scaling,
        include_title=True,
        include_motif_plots=include_motif_plots,
        include_bgc_plot=True,
        include_compound_plots=include_compound_plots,
        motif_reports_base_url="../motif_reports",
    )
    
    output_file = bgc_dir / f"{bgc_id}.html"
    with open(output_file, "w") as f:
        f.write(html_content)
    
    motif_ids = [hit["motif_id"] for hit in motif_hits]
    
    return {
        "href": f"{bgc_id}.html",
        "bgc_id": bgc_id,
        "num_motifs": len(motif_ids),
        "motif_ids": motif_ids
    }


def _generate_single_motif_html(args):
    """Worker function to generate a single motif HTML report."""
    (motif_id, bgcs, global_data, motif_dir, gene_arrow_scaling, 
     include_compound_plots, include_motif_plots) = args
    
    motif_hits = global_data.motif2hits.get(motif_id, [])
    
    html_content = generate_html_for_motif(
        motif_id=motif_id,
        motif_hits=motif_hits,
        bgcs=bgcs,
        domains=global_data.domain_hits,
        domain_colors=global_data.domain_colors,
        motifs=global_data.motifs,
        compounds_info=global_data.compounds,
        scaling=gene_arrow_scaling,
        include_title=True,
        include_motif_plots=include_motif_plots,
        include_compound_plots=include_compound_plots,
        bgc_reports_base_url="../bgc_reports",
    )
    
    output_file = motif_dir / f"{motif_id}.html"
    with open(output_file, "w") as f:
        f.write(html_content)
    
    bgc_ids = [hit["bgc_id"] for hit in motif_hits]
    
    return {
        "href": f"{motif_id}.html",
        "motif_id": motif_id,
        "num_bgcs": len(bgc_ids),
        "bgc_ids": bgc_ids
    }


class SubSketchSession:
    """High-level interface for SubSketch with lazy BGC loading."""

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
        self._bgc_cache: Dict[str, dict] = {}  # Single source of truth for BGCs

    def load(self, n_jobs: int = 1, load_bgcs_upfront: bool = False) -> None:
        """Load all data.
        
        Args:
            n_jobs: Number of parallel processes for loading GenBank files (default: 1)
            load_bgcs_upfront: If True, load all GenBank files immediately; 
                            if False, load on-demand (default: False)
        """
        # Load metadata files
        logger.info("Loading metadata files...")
        domain_hits = io.read_domain_hits(self.domain_hits_file)
        bgc2hits, motif2hits = io.read_detected_motifs(self.motif_hits_file)
        motifs = io.read_motifs(self.motifs_file)
        compounds = io.read_compounds(self.compounds_file) if self.compounds_file else {}

        domain_colors = loaders.load_domain_colors(self.domain_colors_file)
        new_colors = loaders.new_domain_colors(domain_hits, domain_colors)
        if new_colors:
            logger.info(f"  Generating colors for {len(new_colors)} new protein domains")
            domain_colors.update(new_colors)
            io.write_domain_colors(domain_colors, self.domain_colors_file)

        self.data = GlobalData(
            domain_hits=domain_hits,
            bgc2hits=bgc2hits,
            motif2hits=motif2hits,
            motifs=motifs,
            compounds=compounds,
            domain_colors=domain_colors,
        )

        if load_bgcs_upfront:
            # Load all BGCs upfront into cache
            gbk_paths = sorted(self.genbank_dir.glob("*.gbk"))
            logger.info(f"Loading {len(gbk_paths)} GenBank files using {n_jobs} processes")
            
            if n_jobs == 1:
                self._bgc_cache = {}
                for i, gbk_path in enumerate(gbk_paths, 1):
                    if i % 1000 == 0:
                        logger.info(f"  Loaded {i}/{len(gbk_paths)} BGCs...")
                    bgc = loaders.load_bgc(gbk_path)
                    self._bgc_cache[bgc["id"]] = bgc
            else:
                with Pool(processes=n_jobs) as pool:
                    bgc_list = pool.map(loaders.load_bgc, gbk_paths)
                self._bgc_cache = {bgc["id"]: bgc for bgc in bgc_list}
            
            logger.info(f"Loaded {len(self._bgc_cache)} BGCs")
        else:
            logger.info("Lazy loading mode: BGCs will be loaded on-demand as needed")
            self._bgc_cache = {}

    def _get_bgc(self, bgc_id: str) -> dict:
        """Get BGC data, loading from disk if necessary."""
        if bgc_id in self._bgc_cache:
            return self._bgc_cache[bgc_id]
        
        # Load from disk
        gbk_path = self.genbank_dir / f"{bgc_id}.gbk"
        if not gbk_path.exists():
            raise KeyError(f"BGC GenBank file not found: {gbk_path}")
        
        logger.debug(f"Lazy loading BGC: {bgc_id}")
        bgc_data = loaders.load_bgc(gbk_path)
        self._bgc_cache[bgc_id] = bgc_data
        return bgc_data

    def _get_multiple_bgcs(self, bgc_ids: List[str], n_jobs: int = 1) -> Dict[str, dict]:
        """Get multiple BGCs efficiently, loading only what's needed."""
        # Separate cached from uncached
        result = {}
        to_load = []
        
        for bgc_id in bgc_ids:
            if bgc_id in self._bgc_cache:
                result[bgc_id] = self._bgc_cache[bgc_id]
            else:
                to_load.append(bgc_id)
        
        if to_load:
            logger.info(f"Loading {len(to_load)} BGCs on-demand...")
            gbk_paths = [self.genbank_dir / f"{bgc_id}.gbk" for bgc_id in to_load]
            
            if n_jobs == 1:
                for gbk_path in gbk_paths:
                    bgc = loaders.load_bgc(gbk_path)
                    result[bgc["id"]] = bgc
                    self._bgc_cache[bgc["id"]] = bgc
            else:
                with Pool(processes=n_jobs) as pool:
                    bgc_list = pool.map(loaders.load_bgc, gbk_paths)
                for bgc in bgc_list:
                    result[bgc["id"]] = bgc
                    self._bgc_cache[bgc["id"]] = bgc
        
        return result

    def list_genbanks(self) -> List[str]:
        """Return all genbanks available in the BGC input dir."""
        if self.data is None:
            raise RuntimeError("Session not loaded. Call .load() first.")
        
        # Return all available BGC IDs (from files)
        return sorted([p.stem for p in self.genbank_dir.glob("*.gbk")])

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
        
        bgc_data = self._get_bgc(bgc_id)
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
        
        # Load only the BGCs needed for this motif
        bgc_ids = [hit["bgc_id"] for hit in motif_hits]
        bgcs = self._get_multiple_bgcs(bgc_ids, n_jobs=1)

        return generate_html_for_motif(
            motif_id=motif_id,
            motif_hits=motif_hits,
            bgcs=bgcs,
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
        n_jobs: int = 1,
    ) -> None:
        """Generate all reports (BGC, Motif) with a master index page.
        
        Args:
            output_dir: Base directory to write all reports
            gene_arrow_scaling: Scaling factor for gene arrows
            include_compound_plots: Whether to include compound visualizations
            include_motif_plots: Whether to include motif plots
            n_jobs: Number of parallel processes (1 = sequential)
        """
        if self.data is None:
            raise RuntimeError("SubSketch Session not loaded. Call .load() first.")

        all_bgc_ids = self.list_genbanks()
        motif_ids = list(self.data.motifs.keys())

        logger.info(
            f"Generating report for {len(all_bgc_ids)} BGCs and {len(motif_ids)} motifs "
            f"unsing {n_jobs} processes"
            )
        
        bgc_dir = Path(output_dir) / "bgc_reports"
        motif_dir = Path(output_dir) / "motif_reports"
        bgc_dir.mkdir(parents=True, exist_ok=True)
        motif_dir.mkdir(parents=True, exist_ok=True)

        # Check for parallel mode with incomplete loading
        if n_jobs > 1 and len(self._bgc_cache) < len(all_bgc_ids):
            logger.error(
                f"Parallel processing requires all GenBank files loaded. "
                f"Currently: {len(self._bgc_cache)}/{len(all_bgc_ids)} BGCs. "
                f"Use n_jobs=1 or reload with load_bgcs_upfront=True")
            raise RuntimeError(
                f"Cannot use parallel processing (n_jobs={n_jobs}) without "
                f"preloading all BGCs. Either set n_jobs=1 or reload session "
                f"with load_bgcs_upfront=True."
            )
        use_parallel = (n_jobs > 1)

        # ========================================
        # Generate BGC reports
        # ========================================
        if use_parallel:
            bgc_args = [
                (bgc_id, self._bgc_cache[bgc_id], self.data, bgc_dir, 
                gene_arrow_scaling, include_compound_plots, include_motif_plots)
                for bgc_id in all_bgc_ids
            ]

            with Pool(processes=n_jobs) as pool:
                bgc_index_entries = []
                for idx, result in enumerate(pool.imap_unordered(
                    _generate_single_bgc_html, bgc_args, chunksize=50), 1):
                    if idx % 1000 == 0:
                        logger.info(f"BGC reports progress: {idx}/{len(all_bgc_ids)}")
                    bgc_index_entries.append(result)

        else:
            bgc_index_entries = []
            for idx, bgc_id in enumerate(all_bgc_ids, 1):
                if idx % 1000 == 0:
                    logger.info(f"BGC reports progress: {idx}/{len(all_bgc_ids)}")
                
                args = (bgc_id, self._get_bgc(bgc_id), self.data, bgc_dir,
                        gene_arrow_scaling, include_compound_plots, include_motif_plots)
                bgc_index_entries.append(_generate_single_bgc_html(args))
            
        logger.info(f"Generated {len(bgc_index_entries)} BGC reports")


        # ========================================
        # Generate motif reports
        # ========================================
        if use_parallel:
            motif_args = []
            for motif_id in motif_ids:
                motif_hits = self.data.motif2hits.get(motif_id, [])
                bgc_ids_needed = [hit["bgc_id"] for hit in motif_hits]
                bgcs_subset = {bid: self._bgc_cache[bid] for bid in bgc_ids_needed}
                motif_args.append((
                    motif_id, 
                    bgcs_subset, 
                    self.data, 
                    motif_dir,
                    gene_arrow_scaling,
                    include_compound_plots,
                    include_motif_plots
                ))
            
            with Pool(processes=n_jobs) as pool:
                motif_index_entries = []
                for idx, result in enumerate(pool.imap_unordered(
                    _generate_single_motif_html, motif_args, chunksize=50), 1):
                    if idx % 1000 == 0:
                        logger.info(f"Motif reports progress: {idx}/{len(motif_ids)}")
                    motif_index_entries.append(result)

        else:
            motif_index_entries = []
            for idx, motif_id in enumerate(motif_ids, 1):
                if idx % 100 == 0:
                    logger.info(f"Motif reports progress: {idx}/{len(motif_ids)}")
                
                motif_hits = self.data.motif2hits.get(motif_id, [])
                bgc_ids_needed = [hit["bgc_id"] for hit in motif_hits]
                bgcs_subset = self._get_multiple_bgcs(bgc_ids_needed, n_jobs=1)
                
                args = (motif_id, bgcs_subset, self.data, motif_dir,
                        gene_arrow_scaling, include_compound_plots, include_motif_plots)
                motif_index_entries.append(_generate_single_motif_html(args))
            
        logger.info(f"Generated {len(motif_index_entries)} motif reports")

        # ========================================
        # Generate master index
        # ========================================
        logger.info("Generating master index...")
        master_index_html = generate_master_index_html(
            bgc_entries=bgc_index_entries,
            motif_entries=motif_index_entries,
            bgc_reports_url="bgc_reports",
            motif_reports_url="motif_reports"
        )
        
        master_index_path = output_dir / "index.html"
        with open(master_index_path, "w") as f:
            f.write(master_index_html)
        
        logger.info(f"Reports generated successfully. Open {master_index_path} to view.")