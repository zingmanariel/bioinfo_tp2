"""Shared plots: the dot-plot, the path drawn on top of it, and the protein
comparison matrix.

Plots used by a single exercise live in that exercise's own module.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

BLUE = '#4a6fa5'
ORANGE = '#d1852f'
GREEN = '#2e8b57'
GREY = '#b0b0b0'
GRID = '#e5e5e5'

# 0 = mismatch (background), 1 = match (dot)
DOT_CMAP = ListedColormap(['white', BLUE])

# Up to this length the sequence letters are written on the axes.
MAX_TICK_LABELS = 50


def _label_axes(ax, dotplot, s1, s2):
    """Axes: the letters if the sequences are short, the indices otherwise."""
    rows, cols = np.asarray(dotplot).shape

    if s1 is not None and len(s1) == rows and rows <= MAX_TICK_LABELS:
        ax.set_yticks(range(rows))
        ax.set_yticklabels(list(s1), fontsize=7)
    if s2 is not None and len(s2) == cols and cols <= MAX_TICK_LABELS:
        ax.set_xticks(range(cols))
        ax.set_xticklabels(list(s2), fontsize=7)

    ax.set_ylabel('s1' if s1 is None else f's1 ({rows})')
    ax.set_xlabel('s2' if s2 is None else f's2 ({cols})')


def _draw_dotplot(ax, dotplot, s1=None, s2=None, title=None):
    array = np.asarray(dotplot, dtype=float)
    ax.imshow(array, cmap=DOT_CMAP, vmin=0, vmax=1,
              interpolation='nearest', aspect='auto')
    if title:
        ax.set_title(title, fontsize=10)
    _label_axes(ax, array, s1, s2)
    return ax


def _draw_path(ax, path, color=ORANGE, label=None):
    xs = [j for i, j in path]
    ys = [i for i, j in path]
    ax.plot(xs, ys, color=color, linewidth=1.6, zorder=3, label=label)
    ax.scatter([xs[0], xs[-1]], [ys[0], ys[-1]], color=color, s=18, zorder=4)


def plot_dotplot(dotplot, s1=None, s2=None, title='Dot-plot'):
    """Draws the dot-plot (3.1b): s1 on the vertical axis, s2 on the horizontal."""
    fig, ax = plt.subplots(figsize=(7, 7))
    _draw_dotplot(ax, dotplot, s1, s2, title)
    fig.tight_layout()
    plt.show()


def plot_dotplot_with_path(dotplot, path, s1=None, s2=None,
                           title='Dot-plot + path'):
    """The dot-plot with the path drawn on top (3.1e)."""
    fig, ax = plt.subplots(figsize=(7, 7))
    _draw_dotplot(ax, dotplot, s1, s2, title)
    _draw_path(ax, path)
    fig.tight_layout()
    plt.show()


def plot_dotplot_grid(dotplots, title='Filter comparison', paths=None):
    """Several dot-plots in a grid, to compare filters or thresholds (3.1c).

    dotplots: dict {label: matrix}. paths: optional dict {label: path}.
    """
    labels = list(dotplots)
    columns = 2 if len(labels) > 1 else 1
    rows = (len(labels) + columns - 1) // columns

    fig, axes = plt.subplots(rows, columns, figsize=(6 * columns, 5.5 * rows))
    axes = np.atleast_1d(axes).ravel()

    for ax, label in zip(axes, labels):
        _draw_dotplot(ax, dotplots[label], title=label)
        if paths and label in paths:
            _draw_path(ax, paths[label])
    for ax in axes[len(labels):]:
        ax.axis('off')

    fig.suptitle(title)
    fig.tight_layout()
    plt.show()


def plot_identity_sweep(true_identities, series, save_path=None,
                        title='Measured identity vs. true identity'):
    """Line plot comparing algorithms across a sweep of true identity values.

    series: list of (label, values, color) tuples, one line each, `values`
    aligned with `true_identities`. Draws a dashed y=x reference line: a
    method that recovers the true identity exactly would sit right on it.
    """
    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    lo, hi = min(true_identities), max(true_identities)
    ax.plot([lo, hi], [lo, hi], linestyle='--', color=GREY, linewidth=1,
            label='y = x (true identity)', zorder=1)

    for label, values, color in series:
        ax.plot(true_identities, values, marker='o', color=color,
                linewidth=1.8, label=label, zorder=2)

    ax.set_xlabel('true identity')
    ax.set_ylabel('measured identity')
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()


def plot_paths_comparison(matrix, paths, s1=None, s2=None,
                          title='Path comparison', diverging=False,
                          save_path=None):
    """Overlays several paths on the same background matrix, to compare
    algorithms (e.g. the greedy path against the DP path).

    paths: list of (label, path, color) tuples.
    diverging=True draws the background with the diverging BLOSUM-style
    colormap used by plot_comparison_matrix(); diverging=False uses the
    binary match/mismatch colormap used by plot_dotplot().
    save_path=None shows the figure interactively; otherwise it's saved to
    that path (and closed) instead of shown.
    """
    array = np.asarray(matrix, dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 7))

    if diverging:
        limit = max(abs(array.min()), abs(array.max())) or 1.0
        image = ax.imshow(array, cmap='RdBu', vmin=-limit, vmax=limit,
                          interpolation='nearest', aspect='auto')
        fig.colorbar(image, ax=ax, shrink=0.8, label='score')
    else:
        ax.imshow(array, cmap=DOT_CMAP, vmin=0, vmax=1,
                  interpolation='nearest', aspect='auto')

    ax.set_title(title, fontsize=10)
    _label_axes(ax, array, s1, s2)

    for label, path, color in paths:
        _draw_path(ax, path, color=color, label=label)
    ax.legend(loc='upper right', fontsize=8)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
    else:
        plt.show()


def plot_comparison_matrix(comparison, s1=None, s2=None,
                           title='Comparison matrix', path=None):
    """Heatmap of the protein score matrix (3.2a / 3.2d).

    Diverging scale centered on 0: blue for positive-scoring pairs
    (conservative substitutions), red for negative ones.
    """
    array = np.asarray(comparison, dtype=float)
    limit = max(abs(array.min()), abs(array.max())) or 1.0

    fig, ax = plt.subplots(figsize=(7.5, 7))
    image = ax.imshow(array, cmap='RdBu', vmin=-limit, vmax=limit,
                      interpolation='nearest', aspect='auto')
    ax.set_title(title, fontsize=10)
    _label_axes(ax, array, s1, s2)
    if path:
        _draw_path(ax, path)
    fig.colorbar(image, ax=ax, shrink=0.8, label='score')
    fig.tight_layout()
    plt.show()
