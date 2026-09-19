"""3.2e) Modifique su algoritmo de busqueda de caminos para que trabaje sobre el
nuevo dot-plot (el de valores reales de 3.2d).

Punto opcional del TP.
"""

from alignment import align_from_path, print_alignment
from dotplot import (build_comparison_matrix, comparison_to_dotplot,
                     comparison_to_scored_dotplot)
from paths import find_path, path_score, step_counts
from plots import plot_comparison_matrix
from sequences import load_protein_pair
from substitution import load_matrix

WINDOW = 7
THRESHOLD = 12


def run():
    s1, s2 = load_protein_pair()
    matrix = load_matrix()
    comparison = build_comparison_matrix(s1, s2, matrix)

    scored = comparison_to_scored_dotplot(comparison, WINDOW, THRESHOLD)
    binary = comparison_to_dotplot(comparison, WINDOW, THRESHOLD)

    scored_path = find_path(scored)
    binary_path = find_path(binary)

    # Both paths are measured with the same yardstick: the BLOSUM score they
    # add up over the comparison matrix (the score over the dot-plot itself
    # isn't comparable, since one is in "cells" and the other in score units).
    for label, path in (('over real values', scored_path),
                        ('over 1s and 0s', binary_path)):
        counts = step_counts(path)
        print(f'path over the {label} dot-plot: {len(path)} cells, '
              f'D={counts["D"]} V={counts["V"]} H={counts["H"]}, '
              f'BLOSUM score of the path {path_score(comparison, path):.0f}')

    row1, row2 = align_from_path(s1, s2, scored_path)
    print('\nalignment from the path over real values:')
    print_alignment(row1, row2)

    plot_comparison_matrix(scored, s1, s2, path=scored_path,
                           title='Path over the dot-plot with real values')
    return scored_path
