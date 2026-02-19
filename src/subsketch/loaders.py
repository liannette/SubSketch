import csv
import sys
from pathlib import Path
from random import uniform
from colorsys import hsv_to_rgb
from subsketch.classes import CodingSequence
from subsketch.io import read_bgc_gbk
from subsketch.utils import extract_domain_base_name


def _extract_cds_features(bgc_id, bgc_record):
    multiple_cds = dict()
    for f in bgc_record.features:
        if f.type == "CDS":
            orf_num = len(multiple_cds) + 1
            multiple_cds[f"{bgc_id}_{orf_num}"] = CodingSequence.from_feature(bgc_id, orf_num, f)
    return multiple_cds


def _extract_mibig_description(bgc_record):
    for f in bgc_record.features:
        if f.type == "misc_feature":
            return f.qualifiers["note"][0]


def _extract_mibig_organism(bgc_record):
    for f in bgc_record.features:
        if f.type == "source":
            return f.qualifiers["organism"][0]


def load_bgc(bgc_gbk_filepath):
    bgc_record = read_bgc_gbk(bgc_gbk_filepath)
    bgc_id = bgc_gbk_filepath.stem
    cds_features = _extract_cds_features(bgc_id, bgc_record)
    bgc_length = len(bgc_record)
    return {
        "id": bgc_id,
        "cds_features": cds_features,
        "length": bgc_length,
        "record": bgc_record,
    }

def load_mibig_bgc(mibig_gbk_filepath):
    bgc_record = read_bgc_gbk(mibig_gbk_filepath)
    bgc_id = bgc_record.annotations.get("locus", "unknown_bgc_id")
    cds_features = _extract_cds_features(bgc_id, bgc_record)
    bgc_length = len(bgc_record)
    description = _extract_mibig_description(bgc_record)
    organism = _extract_mibig_organism(bgc_record)
    return {
        "id": bgc_id,
        "cds_features": cds_features,
        "length": bgc_length,
        "description": description,
        "organism": organism,
        "record": bgc_record,
    }


def load_domain_colors(domains_color_file: str | Path) -> dict:
    """
    Reads and parses a color domains file.

    This function reads a tab-delimited file containing domain names and their
    associated RGB color values. Each line in the file should contain a domain
    name followed by a comma-separated RGB value.

    Args:
        domains_color_file (str): Path to the domains color file.

    Returns:
        dict: A dictionary mapping domain accessions (str) to RGB color values (list of 3 integers).

    Raises:
        SystemExit: If the specified file does not exist.

    Example format of domains_color_file:
        Domain1    255,0,0
        Domain2    0,255,0
        Domain3    0,0,255
    """
    if not Path(domains_color_file).is_file():
        sys.exit(f"Error while loading domain colors. File was not found: {domains_color_file}")

    domain_colors = dict()
    with open(domains_color_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            domain_accession = row[0]
            domain_accession = extract_domain_base_name(domain_accession)
            rgb = [int(val) for val in row[1].split(",")]
            domain_colors[domain_accession] = rgb
    return domain_colors


def new_domain_colors(domain_hits, existing_domain_colors) -> None:
    """
    Return colors ONLY for domains not yet in existing_domain_colors.
    """
    all_domains = {
        extract_domain_base_name(domain_hit["accession"])
        for bgc_id, orf_hits in domain_hits.items()
        for orf_number, domain_hits_list in orf_hits.items()
        for domain_hit in domain_hits_list
    }

    new_colors = dict()
    for acc in all_domains:
        if acc not in existing_domain_colors:
            h = uniform(0, 1)  # Random hue
            s = uniform(0.5, 0.75)  # Saturation (less)
            v = uniform(0.7, 0.9)  # Value (darker)
            r, g, b = tuple(round(c * 255) for c in hsv_to_rgb(h, s, v))
            new_colors[acc] = (r, g, b)

    return new_colors
    
