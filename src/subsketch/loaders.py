
from subsketch.classes import CodingSequence
from subsketch.io import read_bgc_gbk

def _extract_cds_features(bgc_id, bgc_record):
    multiple_cds = dict()
    for f in bgc_record.features:
        if f.type == "CDS":
            orf_num = len(multiple_cds) + 1
            multiple_cds[f"{bgc_id}_{orf_num}"] = CodingSequence(bgc_id, orf_num, f)
    return multiple_cds


def _extract_mibig_description(bgc_record):
    for f in bgc_record.features:
        if f.type == "misc_feature":
            return f.qualifiers["note"][0]


def _extract_mibig_organism(bgc_record):
    for f in bgc_record.features:
        if f.type == "source":
            return f.qualifiers["organism"][0]


def load_bgc_data(bgc_gbk_filepath):
    bgc_record = read_bgc_gbk(bgc_gbk_filepath)
    bgc_id = bgc_record.annotations.get("locus", "unknown_bgc_id")
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