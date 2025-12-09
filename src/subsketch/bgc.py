# Authors: Annette Lien (2025), Joris Louwen (2019), Jorge Navarro (2016), Peter Cimermancic (2010)

from math import atan2, pi, sin
from colorsys import hsv_to_rgb, rgb_to_hsv

from subsketch.config import (
    internal_domain_margin,
    domain_contour_thickness,
    gene_contour_thickness,
    stripe_thickness,
)

def _get_gene_coordinates(X, Y, L, l, H, h, strand):
    if strand == "+":
        if L < l:
            # squeeze arrow if length shorter than head length
            A = [X, Y - h]
            B = [X + L, Y + H / 2]
            C = [X, Y + H + h]
            points = [A, B, C]
        else:
            A = [X, Y]
            B = [X + L - l, Y]
            C = [X + L - l, Y - h]
            D = [X + L, Y + H / 2]
            E = [X + L - l, Y + H + h]
            F = [X + L - l, Y + H]
            G = [X, Y + H]
            points = [A, B, C, D, E, F, G]

    elif strand == "-":
        if L < l:
            # squeeze arrow if length shorter than head length
            A = [X, Y + H / 2]
            B = [X + L, Y - h]
            C = [X + L, Y + H + h]
            points = [A, B, C]
        else:
            A = [X + L, Y]
            B = [X + l, Y]
            C = [X + l, Y - h]
            D = [X, Y + H / 2]
            E = [X + l, Y + H + h]
            F = [X + l, Y + H]
            G = [X + L, Y + H]
            points = [A, B, C, D, E, F, G]
    return points


def _get_arrow_head_location(L, l, strand):
    if strand == "+":
        head_end = L
        # no tail
        if L < l:
            head_start = 0
        # tail
        else:
            head_start = L - l  # relative to start of gene, not absolute coords.
    elif strand == "-":
        head_start = 0
        # no tail
        if L < l:
            head_end = L
        # tail
        else:
            head_end = l

    return head_start, head_end


def _get_domain_coordinates(
    dX, dL, dH, X, Y, L, H, h, strand, head_start, head_end, head_length
):
    """
    Get coordinates for a domain on a gene arrow.

    Domains on the tip of the arrow should not have corners sticking out
    (which would happen if we just drew a rectangle).

    Args:
        dX: domain start relative to gene start
        dL: domain length
        dH: domain height
        X: gene start x coordinate
        Y: gene start y coordinate
        L: gene length
        H: gene height
        h: gene head edge height
        strand: gene strand ("+" or "-")
        head_start: position where arrow head starts (relative to gene start)
        head_end: position where arrow head ends (relative to gene start)
        head_length: length of the arrow head

    Returns:
        points: list of [x,y] coordinates for the domain polygon
    """
    points = []
    if strand == "+":
        # calculate how far from head_start we would crash with the slope
        # (the horizontal guide at y=Y+internal_domain_margin)
        # Using similar triangles:
        collision_x = head_length * (h + internal_domain_margin)
        collision_x /= h + H / 2.0
        collision_x = round(collision_x)

        # either option for x_margin_offset work
        # m = -float(h + H/2)/(head_length) #slope of right line
        # x_margin_offset = (internal_domain_margin*sqrt(1+m*m))/m
        # x_margin_offset = -(x_margin_offset)
        x_margin_offset = round(
            internal_domain_margin / sin(pi - atan2(h + H / 2.0, -head_length))
        )

        # no collision -> nice, blocky domains
        if (dX + dL) < head_start + collision_x - x_margin_offset:
            points.extend(
                [
                    [X + dX, Y + internal_domain_margin],
                    [X + dX + dL, Y + internal_domain_margin],
                    [X + dX + dL, Y + internal_domain_margin + dH],
                    [X + dX, Y + internal_domain_margin + dH],
                ]
            )
        # collision -> draw a polygon
        else:
            points = []

            # handle the left part of domain (tail)
            if dX < head_start + collision_x - x_margin_offset:
                # arrow with tail: add points A and B
                points.append([X + dX, Y + internal_domain_margin])
                points.append(
                    [
                        X + head_start + collision_x - x_margin_offset,
                        Y + internal_domain_margin,
                    ]
                )
            else:
                # arrow without tail: add point A'
                start_y_offset = int(
                    (h + H / 2) * (L - x_margin_offset - dX) / head_length
                )
                points.append([X + dX, int(Y + H / 2 - start_y_offset)])

            # handle the right part of domain (arrow head)
            if dX + dL >= head_end - x_margin_offset:  # could happen with scaling
                # right part is a triangle
                points.append([X + head_end - x_margin_offset, int(Y + H / 2)])
            else:
                # right part is a cut triangle
                end_y_offset = (2 * h + H) * (L - x_margin_offset - dX - dL)
                end_y_offset /= 2 * head_length
                end_y_offset = int(end_y_offset)
                points.append([X + dX + dL, int(Y + H / 2 - end_y_offset)])
                points.append([X + dX + dL, int(Y + H / 2 + end_y_offset)])

            # handle lower part
            if dX < head_start + collision_x - x_margin_offset:
                # arrow with tail: add points E and F
                points.append(
                    [
                        X + head_start + collision_x - x_margin_offset,
                        Y + H - internal_domain_margin,
                    ]
                )
                points.append([X + dX, Y + H - internal_domain_margin])
            else:
                # # arrow without tail: add point F'
                points.append([X + dX, int(Y + H / 2 + start_y_offset)])

    # now check other direction (strand == "-")
    elif strand == "-":
        # calculate how far from head_start we would crash with the slope
        # (the horizontal guide at y=Y+internal_domain_margin)
        # Using similar triangles:
        collision_x = head_length * ((H / 2) - internal_domain_margin)
        collision_x /= h + H / 2.0
        collision_x = round(collision_x)

        x_margin_offset = round(
            internal_domain_margin / sin(atan2(h + H / 2.0, head_length))
        )

        # no collision -> nice, blocky domains
        if dX > collision_x + x_margin_offset:
            points.extend(
                [
                    [X + dX, Y + internal_domain_margin],
                    [X + dX + dL, Y + internal_domain_margin],
                    [X + dX + dL, Y + internal_domain_margin + dH],
                    [X + dX, Y + internal_domain_margin + dH],
                ]
            )
        # collision -> draw a polygon
        else:
            # handle left part of domain (head)
            if dX < x_margin_offset:
                # regular triangle
                points.append([X + x_margin_offset, Y + H / 2])
            else:
                # cut triangle
                start_y_offset = round(
                    (h + H / 2) * (dX - x_margin_offset) / head_length
                )
                points.append([X + dX, Y + H / 2 + start_y_offset])
                points.append([X + dX, Y + H / 2 - start_y_offset])

            # handle middle/end
            if dX + dL < collision_x + x_margin_offset:
                # no tail
                if head_length != 0:
                    end_y_offset = round(
                        (h + H / 2) * (dX + dL - x_margin_offset) / head_length
                    )
                else:
                    end_y_offset = 0
                points.append([X + dX + dL, Y + H / 2 - end_y_offset])
                points.append([X + dX + dL, Y + H / 2 + end_y_offset])
            else:
                # tail
                points.append(
                    [X + collision_x + x_margin_offset, Y + internal_domain_margin]
                )
                points.append([X + dX + dL, Y + internal_domain_margin])
                points.append([X + dX + dL, Y + internal_domain_margin + dH])
                points.append(
                    [
                        X + collision_x + x_margin_offset,
                        Y + internal_domain_margin + dH,
                    ]
                )
    return points


def _get_domain_stroke_rgb(domain_fill_rgb):
    # contour color is a bit darker. We go to h,s,v space for that
    h_, s, v = rgb_to_hsv(
        float(domain_fill_rgb[0]) / 255.0,
        float(domain_fill_rgb[1]) / 255.0,
        float(domain_fill_rgb[2]) / 255.0,
    )
    domain_stroke_rgb = tuple(int(c * 255) for c in hsv_to_rgb(h_, s, 0.8 * v))
    return domain_stroke_rgb


# --- Draw arrow for gene
def draw_arrow(
    additional_tabs,
    X,
    Y,
    L,
    l,  # noqa: E741, E741
    H,
    h,
    strand,
    color,
    color_contour,
    cds_tag,
    domain_list,
    domain_colors,
    scaling,
    arrow_opacity=1.0,
    domain_opacity=0.75,
):
    """
    SVG code for arrow:
        - (X,Y) ... upper left (+) or right (-) corner of the arrow
        - L ... arrow length
        - H ... arrow height
        - strand
        - h ... arrow head edge width
        - l ... arrow head length
        - color
        - strand
    the edges are ABCDEFG starting from (X,Y)
    domain_list: list of elements to draw domains
    """
    if strand not in ["+", "-"]:
        return ""

    arrow_head_start, arrow_head_end = _get_arrow_head_location(L, l, strand)
    arrow_head_length = arrow_head_end - arrow_head_start

    if arrow_head_length == 0:
        return ""

    svg_str = f"{additional_tabs}<g>\n"

    gene_points = _get_gene_coordinates(X, Y, L, l, H, h, strand)

    svg_str += f"{additional_tabs}\t<title>{cds_tag}</title>\n"
    svg_str += (
        f'{additional_tabs}\t<polygon class="{cds_tag}" '
        f'points="{" ".join(f"{point[0]},{point[1]}" for point in gene_points)}" '
        f'fill="rgb({",".join([str(val) for val in color])})" fill-opacity="1.0" '
        f'stroke="rgb({",".join([str(val) for val in color_contour])})" '
        f'stroke-width="{gene_contour_thickness}" opacity="{arrow_opacity}"/>\n'
    )

    domain_height = int(H - 2 * internal_domain_margin)
    for domain in domain_list:
        domain_accession = domain["accession"]
        domain_start = int(domain["start"] / scaling)
        domain_width = int(domain["width"] / scaling)
        domain_fill_rgb = domain_colors[domain_accession]
        domain_stroke_rgb = _get_domain_stroke_rgb(domain_fill_rgb)
        domain_points = _get_domain_coordinates(
            domain_start,
            domain_width,
            domain_height,
            X,
            Y,
            L,
            H,
            h,
            strand,
            arrow_head_start,
            arrow_head_end,
            arrow_head_length,
        )
        domain_points_str = " ".join(f"{p[0]},{p[1]}" for p in domain_points)

        svg_str += f"{additional_tabs}\t<g>\n"
        svg_str += f'{additional_tabs}\t\t<title>{domain["accession"]}"</title>\n'
        svg_str += (
            f'{additional_tabs}\t\t<polygon class="{domain["accession"]}" '
            f'points="{domain_points_str}" stroke-linejoin="round" '
            f'width="{domain_width}" height="{domain_height}" '
            f'fill="rgb({",".join([str(val) for val in domain_fill_rgb])})" '
            f'stroke="rgb({",".join([str(val) for val in domain_stroke_rgb])})" '
            f'stroke-width="{domain_contour_thickness}" opacity="{domain_opacity}" />\n'
        )
        svg_str += f"{additional_tabs}\t</g>\n"

    # end of gene arrow
    svg_str += f"{additional_tabs}</g>\n"

    return svg_str


def draw_line(X, Y, L):
    """
    Draw a line below genes
    """
    return (
        f'<line x1="{X}" y1="{Y}" x2="{X + L}" y2="{Y}" '
        f'style="stroke:rgb(70,70,70); stroke-width:{stripe_thickness} "/>\n'
    )


def _get_tokenized_gene(domains):
    return ";".join([domain["accession"] for domain in domains])


def draw_bgc(
    bgc_data,
    highlighted_cds=set(),
    bgc_domain_hits=dict(),
    domain_colors=dict(),
    color_genes=True,
    color_domains=True,
    H=30,
    h=5,
    l=12,  # noqa: E741
    mX=10,
    mY=10,
    scaling=30,
    additional_tabs="\t",
):
    """
    Renders a BGC (Biosynthetic Gene Cluster) or detected motif as SVG visualization.

    This function creates an SVG representation of a full biosynthetic gene cluster of
    a specific motif. It draws gene arrows with domain information and includes
    relevant annotations from the GenBank file.

    Args:
        bgc_length (int): Length of the BGC in base pairs.
        cds_features (dict): Dictionary mapping cds_id to CodingSequence objects.
            Format: {"{bgc_id}_{cds_num}": CodingSequence}.
        bgc_domain_hits (dict): Dictionary mapping orf_num to domain hits.
        domain_colors (dict): Dictionary mapping domain accessions to RGB color tuples.
        motif_hit (dict, optional): Information about a detected motif. If provided, only genes
            matching the motif will be drawn. Dictionary should contain keys:
            - motif_id: Identifier for the motif
            - n_matches: Number of matches found
            - threshold: Detection threshold
            - score: Motif score
            - genes: List of genes belonging to the motif
        H (int, optional): Height of gene arrows in SVG units. Defaults to 30.
        h (int, optional): Height of gene arrow head edge in SVG units. Defaults to 5.
        l (int, optional): Maximum length of gene arrow head in SVG units. Defaults to 12.
        mX (int, optional): Horizontal padding in SVG units. Defaults to 10.
        mY (int, optional): Vertical padding in SVG units. Defaults to 10.
        scaling (int, optional): Horizontal scaling factor to convert sequence positions to SVG coordinates.
            Larger values create a more compact visualization. Defaults to 30.

    Returns:
        str: SVG markup as a string that can be embedded in HTML or saved to a file.
    """
    # Calculate SVG dimensions
    svg_width = bgc_data["length"] / scaling + 2 * mX
    svg_height = 2 * h + H + 2 * mY
    svg_text = f'<svg width="{svg_width}" height="{svg_height}">\n'

    # draw a line that corresponds to cluster size
    line = draw_line(mX, mY + h + H / 2, bgc_data["length"] / scaling)
    svg_text += f"{additional_tabs}{line}"

    # gene colors
    color_fill_as = {
        "biosynthetic": (129, 14, 21),
        "biosynthetic-additional": (241, 109, 117),
        "transport": (241, 109, 117),
        "regulatory": (46, 139, 87),
        "other": (128, 128, 128),
    }

    # draw genes, one by one
    for cds in bgc_data["cds_features"].values():

        # standard colors (black contour, white fill)
        color_contour = (0, 0, 0)
        color_fill = (255, 255, 255)
        # color by antiSMASH gene type annotation
        if color_genes:
            color_fill = color_fill_as.get(cds.as_type, (255, 255, 255))

        # Set gene arrow opacity and protein domains in gene
        if highlighted_cds and cds.orf_num not in highlighted_cds:
            arrow_opacity = 0.2
            domain_opacity = 0.1
        else:
            arrow_opacity = 1.0
            domain_opacity = 0.75
        
        domains = bgc_domain_hits.get(cds.orf_num, []) if color_domains else []

        # define arrow's start and end
        # http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html#start
        gene_start = cds.start / scaling
        gene_end = cds.end / scaling
        gene_length = gene_end - gene_start

        arrow = draw_arrow(
            additional_tabs=additional_tabs,
            X=gene_start + mX,
            Y=mY + h,
            L=gene_length,
            l=l,
            H=H,
            h=h,
            strand=cds.strand,
            color=color_fill,
            color_contour=color_contour,
            cds_tag=cds.tag,
            domain_list=domains,
            domain_colors=domain_colors,
            scaling=scaling,
            arrow_opacity=arrow_opacity,
            domain_opacity=domain_opacity,
        )
        if arrow == "":
            print(f"  (ArrowerSVG) Warning: something went wrong with {cds.bgc_id}_orf{cds.orf_num}")
        svg_text += arrow

    # Close the SVG tag
    svg_text += "</svg>\n"

    return svg_text


def draw_subcluster_hit(
    bgc_data,
    bgc_domain_hits,
    domain_colors,
    motif_hit,
    H=30,
    h=5,
    l=12,  # noqa: E741
    mX=10,
    mY=10,
    scaling=30,
):
    # Filter cds that are not part of the detected motif
    highlighted_cds = set()
    for cds in bgc_data["cds_features"].values():
        domains = bgc_domain_hits.get(cds.orf_num, [])
        if _get_tokenized_gene(domains) in motif_hit["genes"]:
            highlighted_cds.add(cds.orf_num)

    svg_text = draw_bgc(
        bgc_data=bgc_data,
        highlighted_cds=highlighted_cds,
        bgc_domain_hits=bgc_domain_hits,
        domain_colors=domain_colors,
        color_genes=False,
        color_domains=True,
        H=H,
        h=h,
        l=l,
        mX=mX,
        mY=mY,
        scaling=scaling,
    )
    return svg_text


def add_subcluster_title(svg_text, motif_hit):
    """Wrap SVG content with a title in HTML format.

    Args:
        svg_text (str): The SVG content as a string.
        motif_hit (dict): The motif hit information containing details for the title.

    Returns:
        str: HTML string containing the title and SVG content.
    """
    title = f"Subcluster Motif #{motif_hit['motif_id']}"
    subtitle = f"Threshold: {motif_hit['threshold']}\nScore: {motif_hit['score']}"

    html_content = f"""
    <div style="margin: 1.5em 0 1em 0;">
        <h3 style="margin: 0 0 0.3em 0; line-height: 1.2;">
            {title}
            <br><small style="font-size: 0.8em; color: #666; font-weight: 400;">
                {subtitle}
            </small>
        </h3>
    </div>
    <div>{svg_text}</div>
    """
    return html_content


def draw_annotated_subcluster(
    annotated_subcluster,
    bgc_data,
    bgc_domain_hits=dict(),
    domain_colors=dict(),
    color_genes=True,
    color_domains=False,
    H=30,
    h=5,
    l=12,  # noqa: E741
    mX=10,
    mY=10,
    scaling=30,
):
    # Filter cds that are not part of the annotated motif
    highlighted_cds = set()
    for cds in bgc_data["cds_features"].values():
        protein_id = cds.seqio_feature.qualifiers["protein_id"][0]
        if protein_id in annotated_subcluster["protein_ids"]:
            highlighted_cds.add(cds.orf_num)

    svg_text = draw_bgc(
        bgc_data=bgc_data,
        highlighted_cds=highlighted_cds,
        bgc_domain_hits=bgc_domain_hits,
        domain_colors=domain_colors,
        color_genes=color_genes,
        color_domains=color_domains,
        H=H,
        h=h,
        l=l,
        mX=mX,
        mY=mY,
        scaling=scaling,
    )
    return svg_text