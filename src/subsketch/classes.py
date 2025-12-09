import re
import numpy as np
from typing import TYPE_CHECKING, Optional, List
from dataclasses import dataclass

if TYPE_CHECKING:
    from Bio.SeqFeature import SeqFeature

@dataclass
class CodingSequence:
    bgc_id: str
    orf_num: int
    as_type: str
    protein_id: Optional[str]
    locus_tag: Optional[str]
    gene_name: Optional[str]
    gene_product: Optional[str]
    start: int
    end: int
    strand: str
    seqio_feature: "SeqFeature"

    @classmethod
    def from_feature(cls, bgc_id: str, orf_num: int, seqio_feature: "SeqFeature") -> "CodingSequence":
        """Create from Biopython SeqFeature."""
        qualifiers = seqio_feature.qualifiers
        
        return cls(
            bgc_id=bgc_id,
            orf_num=orf_num,
            as_type=qualifiers.get("gene_kind", ["other"])[0],
            protein_id=qualifiers.get("protein_id", [None])[0],
            locus_tag=qualifiers.get("locus_tag", [None])[0],
            gene_name=qualifiers.get("gene", [None])[0],
            gene_product=qualifiers.get("product", [None])[0],
            start=int(seqio_feature.location.start),
            end=int(seqio_feature.location.end),
            strand=cls._numerical_strand_to_string(seqio_feature.location.strand),
            seqio_feature=seqio_feature
        )
    
    @staticmethod
    def _numerical_strand_to_string(strand: int) -> str:
        """Convert numerical strand annotation to +/- string."""
        if strand == 1:
            return "+"
        elif strand == -1:
            return "-"
        raise ValueError(f"Invalid strand: {strand}")
    
    @property
    def tag(self) -> str:
        """SVG-friendly tag combining all identifiers."""
        parts = [
            f"gene: {self.gene_name}",
            f"locus_tag: {self.locus_tag}",
            f"protein_id: {self.protein_id}",
            f"product: {self.gene_product}"
        ]
        return " | ".join(p for p in parts if p.split(": ", 1)[1] is not None)
    
    @property
    def length(self) -> int:
        """CDS length."""
        return self.end - self.start


@dataclass
class Motif:
    motif_id: str
    n_matches: int
    threshold: float
    tokenized_genes: List[str]
    weight_matrix: np.ndarray

    @classmethod
    def from_lines(cls, lines: List[str]) -> "Motif":
        """Parse motif from 4-line format."""
        # motif_id, n_matches, threshold = re.findall(r"\d+\.\d+|\d+", lines[0])

        # Extract header: motif_id, n_matches, threshold
        header_match = re.match(r"(\S+)\s+(\d+)\s+(\d+\.\d+)", lines[0].strip())
        if not header_match:
            raise ValueError(f"Invalid motif header: {lines[0]}")
        motif_id, n_matches, threshold = header_match.groups()

        # Parse genes and weights
        tokenized_genes = lines[1].strip().split()
        weights_present = [float(x) for x in lines[2].strip().split()]
        weights_absent = [float(x) for x in lines[3].strip().split()]

        # Validate lengths
        if len(tokenized_genes) != len(weights_present) or len(weights_present) != len(weights_absent):
            raise ValueError("Gene/weight lists have different lengths")
        
        weight_matrix = np.array([weights_present, weights_absent]).T

        return cls(
            motif_id=motif_id,
            n_matches=int(n_matches),
            threshold=float(threshold),
            tokenized_genes=tokenized_genes,
            weight_matrix=weight_matrix
        )
    
    @property
    def n_genes(self) -> int:
        """Number of genes in motif."""
        return len(self.tokenized_genes)
    