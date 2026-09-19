"""Dot-plots: the comparison matrix between two sequences and the filters that
make it readable.

Exercise 3.1 a) and c) for DNA, and 3.2 a), b) and d) for proteins.
Everything programmed in this module is a matrix: numpy does the rest.

Convention used throughout the code: s1 goes in the ROWS (vertical axis,
index i) and s2 in the COLUMNS (horizontal axis, index j).
"""

import numpy as np


# ---------------------------------------------------------------------------
# 3.1a) DNA dot-plot: matrix of 1s (match) and 0s (mismatch)
# ---------------------------------------------------------------------------

def build_dotplot(s1, s2):
    """len(s1) x len(s2) matrix with 1 where s1[i] == s2[j] and 0 otherwise."""
    dotplot = np.zeros((len(s1), len(s2)))
    for i, a in enumerate(s1):
        for j, b in enumerate(s2):
            if a == b:
                dotplot[i][j] = 1
    return dotplot


# ---------------------------------------------------------------------------
# 3.1c) Window filters
# ---------------------------------------------------------------------------
#
# The TP's diagram shows independent decisions, which here are parameters of
# the same function:
#
#   window       window size (w x w)
#   threshold    how many matches the window needs to "turn on"
#                (high threshold vs. low threshold)
#   overlapping  overlapping windows (step 1) or adjacent ones (step w)
#   reduce       the dot-plot is REPLACED by the windows' one (a smaller
#                matrix, "filter and reduce") or the windows are only used as
#                a filter and the matrix keeps its size ("filter, no reduce")
#   diagonal     what gets counted inside the window: every match in the
#                w x w block (what the TP's diagram draws) or only the ones
#                on its diagonal (what a classic dot-plot filter does; see
#                the note below)
#
# Note on `diagonal`: counting the whole block, a window sitting on the true
# diagonal has w real matches + noise, and any random DNA window has on
# average w*w/4 background matches. For w=10 that's 10 + 22 against 25:
# almost no contrast, which is why the block filter barely cleans anything up.
# Counting only the window's diagonal compares w matches against w/4, which is
# the only way the filter actually de-clutters the plot. Exercise 3.1c runs
# both variants and shows the numbers.
# ---------------------------------------------------------------------------

def _window_starts(size, window, overlapping):
    """Coordinate where each window starts along an axis of length `size`."""
    if overlapping:
        return list(range(0, max(size - window + 1, 1)))
    return list(range(0, size, window))


def _window_count(block, diagonal):
    """How many matches the window has: the ones on its diagonal or the whole block's."""
    if diagonal:
        return float(np.trace(block))
    return float(block.sum())


def filter_dotplot(dotplot, window, threshold, overlapping=True, reduce=False,
                   diagonal=False):
    """Applies the window filter over a 0/1 dot-plot.

    Returns another matrix: the same size if reduce=False (the original dots
    of the windows that pass the threshold are kept and the rest is turned
    off), or a matrix with one cell per window if reduce=True.
    """
    array = np.asarray(dotplot, dtype=float)
    rows, cols = array.shape
    row_starts = _window_starts(rows, window, overlapping)
    col_starts = _window_starts(cols, window, overlapping)

    if reduce:
        result = np.zeros((len(row_starts), len(col_starts)))
    else:
        result = np.zeros_like(array)

    for r, i in enumerate(row_starts):
        for c, j in enumerate(col_starts):
            block = array[i:i + window, j:j + window]
            if _window_count(block, diagonal) < threshold:
                continue
            if reduce:
                result[r][c] = 1
            elif diagonal:
                # only turns on the window's diagonal, which is what was counted
                offsets = np.arange(min(block.shape))
                result[i + offsets, j + offsets] = array[i + offsets, j + offsets]
            else:
                result[i:i + window, j:j + window] = block
    return result


def filter_variants(dotplot, window, threshold, diagonal=False):
    """The four combinations from the TP's diagram, to compare them at once.

    Returns a dict {label: filtered matrix}.
    """
    variants = {}
    for overlapping in (True, False):
        for reduce in (False, True):
            label = 'overlapping' if overlapping else 'adjacent'
            label += ', reduced' if reduce else ', not reduced'
            variants[label] = filter_dotplot(dotplot, window, threshold,
                                             overlapping=overlapping,
                                             reduce=reduce,
                                             diagonal=diagonal)
    return variants


# ---------------------------------------------------------------------------
# 3.2a) Protein comparison matrix (substitution-matrix scores)
# ---------------------------------------------------------------------------

def build_comparison_matrix(s1, s2, matrix):
    """len(s1) x len(s2) matrix with the substitution score of each pair.

    Same as build_dotplot(), but instead of 0/1 it stores matrix[a][b]:
    numbers, positive for conservative substitutions and negative for the
    rest. `matrix` is what substitution.load_matrix() returns.
    """
    from substitution import score   # local import: avoids an import cycle

    comparison = np.zeros((len(s1), len(s2)))
    for i, a in enumerate(s1):
        for j, b in enumerate(s2):
            comparison[i][j] = score(matrix, a, b)
    return comparison


# ---------------------------------------------------------------------------
# 3.2b) and 3.2d) From the comparison matrix to the dot-plot
# ---------------------------------------------------------------------------

def _window_mask(comparison, window, threshold, overlapping=True, diagonal=True):
    """Boolean matrix: True where the cell's window passes the threshold.

    overlapping=True: one window centered on each cell (step 1), so the
    dot-plot stays aligned with the sequences. overlapping=False: adjacent
    windows, where the whole block turns on or off together.

    diagonal=True sums the scores along the window's diagonal (a run of
    similar residues), which is what makes a protein dot-plot readable;
    diagonal=False sums the whole block.
    """
    array = np.asarray(comparison, dtype=float)
    rows, cols = array.shape
    mask = np.zeros((rows, cols), dtype=bool)
    half = window // 2

    if overlapping:
        for i in range(rows):
            for j in range(cols):
                if diagonal:
                    total = 0.0
                    for k in range(-half, window - half):
                        if 0 <= i + k < rows and 0 <= j + k < cols:
                            total += array[i + k][j + k]
                else:
                    block = array[max(i - half, 0):i + window - half,
                                  max(j - half, 0):j + window - half]
                    total = float(block.sum())
                mask[i][j] = total >= threshold
        return mask

    for i in _window_starts(rows, window, overlapping=False):
        for j in _window_starts(cols, window, overlapping=False):
            block = array[i:i + window, j:j + window]
            total = float(np.trace(block)) if diagonal else float(block.sum())
            if total >= threshold:
                mask[i:i + window, j:j + window] = True
    return mask


def comparison_to_dotplot(comparison, window, threshold, overlapping=True,
                          diagonal=True):
    """Converts the score matrix into a 1s-and-0s dot-plot.

    Parameters the TP asks for: a window-size and a score threshold. Same
    idea as filter_dotplot(), except what accumulates inside the window is a
    score, not a match count.
    """
    return _window_mask(comparison, window, threshold, overlapping,
                        diagonal).astype(float)


def comparison_to_scored_dotplot(comparison, window, threshold, overlapping=True,
                                 diagonal=True):
    """Same as comparison_to_dotplot() but keeping the REAL values.

    Cells whose window doesn't pass the threshold stay at 0; the rest keep
    their score. This is the dot-plot paths.find_path() works with (3.2e).
    """
    mask = _window_mask(comparison, window, threshold, overlapping, diagonal)
    return np.asarray(comparison, dtype=float) * mask


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def density(dotplot):
    """Fraction of cells turned on: useful to see how much a filter
    "de-clutters" the plot before actually looking at it."""
    array = np.asarray(dotplot)
    return float(np.count_nonzero(array)) / array.size


def preview(dotplot, rows=20, cols=60, on='#', off='.'):
    """Top-left corner of the dot-plot as text, to look at it without
    plotting (the TP's "try it with a smaller matrix first")."""
    array = np.asarray(dotplot)
    lines = []
    for row in array[:rows, :cols]:
        lines.append(''.join(on if value else off for value in row))
    return '\n'.join(lines)
