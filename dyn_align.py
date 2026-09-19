"""Your dynamic-programming alignment pseudocode, translated to Python as-is,
with the helpers that were missing.

What was touched:
  # FIXED  wiring changes, without which the file didn't run

The traceback is already fixed: it walks the whole path, reads the direction
of each cell and discounts score_mat's padding.

Convention (the one your code uses, seq1[col] and seq2[row]):
  rows    = seq2   (row index)
  columns = seq1   (col index)
It's the mirror image of the one dotplot.py uses in the rest of the TP.
"""

import os

from alignment import alignment_stats, print_alignment
from dotplot import build_dotplot, build_comparison_matrix, comparison_to_dotplot
from paths import find_path
from plots import GREEN, ORANGE, plot_dotplot, plot_paths_comparison
from sequences import (DISTANT_PROTEIN, data_path, load_fasta, load_pair,
                       load_protein_pair)
from substitution import load_matrix


class DotPlot:
    """The dot-plot of two sequences, with its matrix and its plot."""

    def __init__(self, seq1, seq2, matrix=None):
        self.seq1 = seq1
        self.seq2 = seq2
        # seq2 goes in the rows and seq1 in the columns, as calc_alignment assumes
        if matrix is not None:
            self.matrix = build_comparison_matrix(seq2, seq1, matrix)
        else:
            self.matrix = build_dotplot(seq2, seq1)

    def plot(self, title='Dot-plot'):
        plot_dotplot(self.matrix, self.seq2, self.seq1, title=title)


GAP = -1


def get_best_move(dot_plot, score_mat, row, col, gap):
    """(score, parent) of the best step that reaches (row, col).

    A step's score depends on HOW you get there, so choosing the parent and
    computing the score are the same decision: whichever maximizes parent's
    score + step's score wins.

      DIAG   aligns seq2[row-1] with seq1[col-1] -> whatever the dot-plot says
             (1 if a match, 0 if a mismatch)
      VERT   seq2's letter against a gap        -> GAP
      HORIZ  seq1's letter against a gap        -> GAP

    max() keeps the first of any tied options, and the diagonal comes first,
    so a tie is won by the diagonal step.
    """
    match_score = get_cell_score(dot_plot, row, col)
    options = [
        (score_mat[row - 1][col - 1] + match_score, (row - 1, col - 1)),   # DIAG
        (score_mat[row - 1][col] + gap, (row - 1, col)),                   # VERT
        (score_mat[row][col - 1] + gap, (row, col - 1)),                   # HORIZ
    ]
    return max(options, key=lambda option: option[0])


def get_cell_score(dot_plot, row, col):
    """Dot-plot value for cell (row-1, col-1); 0 if it falls outside."""
    if row - 1 < len(dot_plot.matrix) and col - 1 < len(dot_plot.matrix[0]):
        return dot_plot.matrix[row - 1][col - 1]
    return 0


def get_dir(path, i):
    """Direction of the step that reaches path[i] from its parent.

    path comes in reverse (from the end to the start), so path[i]'s parent is
    path[i + 1]. The last cell's parent is the starting point (0, 0), which
    isn't in path: if the path starts right at the border, that step is a gap.
    """
    row, col = path[i]
    parent_row, parent_col = path[i + 1] if i + 1 < len(path) else (0, 0)
    if parent_row == row - 1 and parent_col == col - 1:
        return "DIAG"
    if parent_row == row - 1:
        return "VERT"
    return "HORIZ"


def calc_alignment(seq1, seq2, dot_plot, gap=GAP):

    rows = len(dot_plot.matrix)
    cols = len(dot_plot.matrix[0])

    score_mat = [[0] * (cols + 1) for _ in range(rows + 1)]
    best_parent_mat = [[None] * (cols + 1) for _ in range(rows + 1)]

    for row in range(len(score_mat)):
        for col in range(len(score_mat[0])):
            if col == 0 or row == 0:
                score_mat[row][col] = 0
                # The border is the initial run of gaps, so it has a parent
                # too: each cell comes from the previous one on the same
                # border. (0, 0) is the path's start and stays None, which is
                # what stops the traceback's while loop.
                if col > 0:
                    best_parent_mat[row][col] = (row, col - 1)
                elif row > 0:
                    best_parent_mat[row][col] = (row - 1, col)
            else:
                # the step's score depends on its direction, so the parent
                # and the score come out of the same comparison
                score_mat[row][col], best_parent_mat[row][col] = get_best_move(dot_plot, score_mat, row, col, gap)  # row-1, col-1, if pos is outside of dotplot score is 0

    #----------
    best_path = []

    row = len(score_mat) - 1
    col = len(score_mat[0]) - 1

    while not(row == 0 and col == 0):
        best_path.append((row, col))
        row, col = best_parent_mat[row][col]

    #-----------
    seq1_aligned = []
    seq2_aligned = []

    # looking from the last position, if we moved then we put the letter from the current possition, if we didnt move we put a gap
    # FIXED: the list is called best_path (it used to say path here)
    # Cell (row, col) of score_mat is cell (row-1, col-1) of the dot-plot,
    # i.e. seq2[row-1] and seq1[col-1].
    for i in range(len(best_path)):
        last_dir = get_dir(best_path, i)
        row, col = best_path[i]

        if last_dir == "DIAG":
            seq1_aligned.insert(0, seq1[col-1])
            seq2_aligned.insert(0, seq2[row-1])

        elif last_dir == "VERT":
            seq1_aligned.insert(0, "-")
            seq2_aligned.insert(0, seq2[row-1])

        else:  # last_dir == "HORIZ"
            seq1_aligned.insert(0, seq1[col-1])
            seq2_aligned.insert(0, "-")

    # best_path is in score_mat coordinates (row-1, col-1) = (seq2 index, seq1
    # index) and runs end-to-start; convert to dot-plot coordinates and put it
    # start-to-end, for anyone who wants to plot it (see save_figures()).
    path = [(row - 1, col - 1) for row, col in reversed(best_path)]

    return ''.join(seq1_aligned), ''.join(seq2_aligned), path


# ---------------------------------------------------------------------------
# Comparison against the TP's proposed algorithm (dotplot + path search) and
# against BLAST: this is where the numbers used in the report come from.
# ---------------------------------------------------------------------------

DNA_GAP = -1
PROTEIN_GAPS = [-1, -4, -8, -11]   # -11 ~ blastp's gap-open with BLOSUM62


def run_dna():
    """Real DNA: human beta-globin vs. mouse beta-globin (same pair 3.1 uses)."""
    s1, s2 = load_pair()
    dot_plot = DotPlot(s1, s2)
    row1, row2, _ = calc_alignment(s1, s2, dot_plot, gap=DNA_GAP)
    print(f'--- DNA: human vs. mouse (gap={DNA_GAP}) ---')
    print_alignment(row1, row2)


def run_proteins():
    """Real protein (BLOSUM62): beta-delta (close) and beta-myoglobin
    (distant), same pair 3.2c uses. First the full alignment at a fixed gap,
    then the gap-cost sweep over both pairs."""
    matrix = load_matrix()
    beta, delta = load_protein_pair()
    myoglobin = load_fasta(data_path(DISTANT_PROTEIN))

    pairs = {
        'beta vs. delta (close)': (beta, delta),
        'beta vs. myoglobin (distant)': (beta, myoglobin),
    }

    for label, (s1, s2) in pairs.items():
        dot_plot = DotPlot(s1, s2, matrix=matrix)
        row1, row2, _ = calc_alignment(s1, s2, dot_plot, gap=PROTEIN_GAPS[-1])
        print(f'\n--- {label} (gap={PROTEIN_GAPS[-1]}) ---')
        print_alignment(row1, row2)

    print('\n--- effect of the gap cost ---')
    print(f'{"pair":<32} {"gap":>5} {"length":>7} {"identity":>10} {"gaps":>6}')
    for label, (s1, s2) in pairs.items():
        dot_plot = DotPlot(s1, s2, matrix=matrix)
        for gap in PROTEIN_GAPS:
            row1, row2, _ = calc_alignment(s1, s2, dot_plot, gap=gap)
            stats = alignment_stats(row1, row2)
            print(f'{label:<32} {gap:>5} {stats["length"]:>7} '
                  f'{stats["identity"]:>10.1%} {stats["gaps"]:>6}')


FIGURES_DIR = 'figures'


def _dp_path_in_dotplot_coords(seq1, seq2, dot_plot, gap):
    """DP path, converted from dot_plot's (seq2, seq1) convention to the
    (seq1, seq2) convention the rest of the TP's plots use."""
    _, _, path = calc_alignment(seq1, seq2, dot_plot, gap=gap)
    return [(j, i) for i, j in path]


def save_figures(out_dir=FIGURES_DIR):
    """Saves the report's comparison figures: the greedy path (paths.find_path)
    and the DP path overlaid on the same background matrix, for each of the
    three pairs the report uses."""
    os.makedirs(out_dir, exist_ok=True)

    # DNA: both paths run on the same raw, unfiltered 0/1 dot-plot.
    s1, s2 = load_pair()
    dotplot = build_dotplot(s1, s2)
    greedy_path = find_path(dotplot)
    dp_path = _dp_path_in_dotplot_coords(s1, s2, DotPlot(s1, s2), DNA_GAP)
    plot_paths_comparison(
        dotplot, [('greedy', greedy_path, ORANGE), ('DP', dp_path, GREEN)],
        s1=s1, s2=s2, title='DNA: human vs. mouse',
        save_path=os.path.join(out_dir, 'dna.png'))

    # Proteins: the background is the raw BLOSUM62 comparison matrix; the
    # greedy path runs on the filtered/binarized dot-plot 3.2c uses (w=7,
    # threshold=12), the DP path runs on the raw matrix directly.
    matrix = load_matrix()
    beta, delta = load_protein_pair()
    myoglobin = load_fasta(data_path(DISTANT_PROTEIN))

    pairs = {
        'protein_close': ('beta vs. delta (close)', beta, delta),
        'protein_distant': ('beta vs. myoglobin (distant)', beta, myoglobin),
    }
    for filename, (title, s1, s2) in pairs.items():
        comparison = build_comparison_matrix(s1, s2, matrix)
        filtered = comparison_to_dotplot(comparison, window=7, threshold=12)
        greedy_path = find_path(filtered)
        dp_path = _dp_path_in_dotplot_coords(
            s1, s2, DotPlot(s1, s2, matrix=matrix), PROTEIN_GAPS[-1])
        plot_paths_comparison(
            comparison, [('greedy', greedy_path, ORANGE), ('DP', dp_path, GREEN)],
            s1=s1, s2=s2, title=title, diverging=True,
            save_path=os.path.join(out_dir, f'{filename}.png'))

    print(f'\nfigures saved to {out_dir}/')


def main():
    run_dna()
    run_proteins()
    save_figures()


if __name__ == "__main__":
    main()
