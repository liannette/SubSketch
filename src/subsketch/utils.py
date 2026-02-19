import re


def extract_domain_base_name(domain: str) -> str:
    """Extract the base name of a domain.

    This function extracts the base name of a protein_domain,
    excluding any '_c' suffix of subPfams.

    Args:
        domain (str): The domain name.

    Returns:
        str: The base name of the domain.
    """
    return re.sub(r"_c\d+$", "", domain)