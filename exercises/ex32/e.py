"""3.2e) Modifique su algoritmo de busqueda de caminos para que trabaje sobre el
nuevo dot-plot (el de valores reales de 3.2d).

Punto opcional del TP.
"""

from alignment import align_from_path, print_alignment
from dotplot import (build_comparison_matrix, comparison_to_dotplot,
                     comparison_to_scored_dotplot)
from paths import find_path, find_path_scored, path_score, step_counts
from plots import plot_comparison_matrix
from sequences import load_protein_pair
from substitution import load_matrix

WINDOW = 7
THRESHOLD = 12
GAP_PENALTY = -4.0


def run():
    s1, s2 = load_protein_pair()
    matrix = load_matrix()
    comparison = build_comparison_matrix(s1, s2, matrix)

    scored = comparison_to_scored_dotplot(comparison, WINDOW, THRESHOLD)
    binary = comparison_to_dotplot(comparison, WINDOW, THRESHOLD)

    scored_path = find_path_scored(scored, gap_penalty=GAP_PENALTY)
    binary_path = find_path(binary)

    # Los dos caminos se miden con la misma vara: el puntaje BLOSUM que suman
    # sobre la matriz de comparacion (el puntaje sobre el propio dot-plot no es
    # comparable, porque uno esta en unidades de "casillas" y el otro en
    # unidades de puntaje).
    for label, path in (('con valores reales', scored_path),
                        ('de 1s y 0s', binary_path)):
        counts = step_counts(path)
        print(f'camino sobre el dot-plot {label}: {len(path)} casillas, '
              f'D={counts["D"]} V={counts["V"]} H={counts["H"]}, '
              f'puntaje BLOSUM del camino {path_score(comparison, path):.0f}')

    row1, row2 = align_from_path(s1, s2, scored_path)
    print(f'\nalineamiento con el camino sobre valores reales '
          f'(gap = {GAP_PENALTY:.0f}):')
    print_alignment(row1, row2)

    plot_comparison_matrix(scored, s1, s2, path=scored_path,
                           title='Camino sobre el dot-plot con valores reales')
    return scored_path
