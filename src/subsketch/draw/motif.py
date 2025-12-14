import io
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def plot_subcluster_motif(motif: dict, hit: dict = None, add_title: bool = False, scaling: float = 1.0):
    """
    Horizontal visualization of subcluster motif hit with stacked bar chart.
    """
    # Extract motif data
    tokens = motif.tokenized_genes
    weight_matrix = np.array(motif.weight_matrix)
    threshold = motif.threshold
    motif_id = motif.motif_id
    pos_widths = weight_matrix[:, 0]   # rightward
    neg_widths = weight_matrix[:, 1]   # leftward

    # Replace inf for negative bars
    INF_REPLACEMENT = 0
    inf_mask = np.isinf(neg_widths)
    neg_widths[inf_mask] = INF_REPLACEMENT

    # Create horizontal plot
    plot_height = len(tokens) * 0.3 * scaling
    plot_width = 6 * scaling
    fig, ax = plt.subplots(figsize=(plot_width, plot_height))

    # Bar heights
    bar_height = 0.8

    # Y positions
    y = np.arange(len(tokens))

    colors = {
        'positive': '#2E8B57',  # Teal
        'negative': '#D2691E',  # Orange
        'muted_pos': '#8AC4A1', # Lighter teal
        'muted_neg': '#E8C8A0'  # Lighter orange
    }


    # Plot positive weights (rightward)
    ax.barh(y, pos_widths, bar_height,
        label='Positive weights', color=colors['muted_pos'], 
        alpha=0.7, linewidth=0.5)
    # Plot negative weights (leftward)
    ax.barh(y, neg_widths, bar_height,
        label='Negative weights', color=colors['muted_neg'], 
        alpha=0.7, linewidth=0.5)

    # Highlight inf negative weights with star
    max_pos = max(weight_matrix[:, 0])
    for i in np.where(inf_mask)[0]:
        ax.text(INF_REPLACEMENT - 0.1 * max_pos, y[i], '★', 
                ha='left', va='center', color=colors["muted_neg"], fontsize=14)
        
    # Highlight matching genes
    if hit is not None:
        hit_genes = set(hit['genes'])
        score = hit['score']
        bgc_id = hit['bgc_id']
        for i, token in enumerate(tokens):
            if token in hit_genes:
                # Navy highlight on positive bar (right side)
                rect = Rectangle((0, y[i] - bar_height/2),
                                 pos_widths[i], bar_height,
                                 facecolor=colors['positive'],
                                 linewidth=0.5, linestyle='-', zorder=10,
                                 alpha=0.8)
            else:
                # Darkred highlight on negative bar (left side)
                rect = Rectangle((neg_widths[i], y[i] - bar_height/2),
                                 abs(neg_widths[i]), bar_height,
                                 facecolor=colors['negative'],
                                 linewidth=0.5, linestyle='-', zorder=10,
                                 alpha=0.8)
            ax.add_patch(rect)
        
    # Title, optional
    if add_title:
        title = f'Subcluster Motif #{motif_id}\n'
        if hit:
            title += f'{bgc_id} - Score: {float(score):.3f} (Threshold: {threshold}) '
        plt.title(title, fontsize=10, fontweight='bold', pad=20)
        plt.subplots_adjust(top=0.88)  # Extra space for titles

    # Y-axis (genes)
    ax.set_yticks(y)
    ax.invert_yaxis()
    # Y-axis labels, bold those in hit_genes
    ytick_texts = ax.set_yticklabels(tokens, fontsize=10)
    if hit is not None:
        hit_genes = set(hit['genes'])
        for text_obj, token in zip(ytick_texts, tokens):
            if token in hit_genes:
                text_obj.set_fontweight('bold')
    
    # X-axis (weights)
    ax.axvline(0, color='black', linewidth=1.2, zorder=5)
    all_widths = np.concatenate([pos_widths, neg_widths])
    ax.set_xlim(left=min(all_widths + [-0.5])*1.2, right=max(pos_widths)*1.2)
    ax.tick_params(axis='x', labelbottom=False)
    ax.set_xlabel('')  

    # Despine
    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")
    ax.set_frame_on(False)

    # return plot as SVG string
    svg_output = io.StringIO()
    plt.savefig(svg_output, format='svg', bbox_inches='tight', dpi=300)
    plt.close(fig)  # Clean up figure
    return svg_output.getvalue()



