"""3.2d) Realice un programa que a partir de la matriz de comparacion y
utilizando filtros obtenga un dot-plot, pero no de 0s y 1s sino que contenga los
valores REALES de la matriz (el filtro convierte en 0 las celdas que no superan
el umbral).

Punto opcional del TP ("si tiene tiempo y ganas").
"""

import numpy as np

from dotplot import (build_comparison_matrix, comparison_to_dotplot,
                     comparison_to_scored_dotplot, density)
from plots import plot_comparison_matrix
from sequences import load_protein_pair
from substitution import load_matrix

WINDOW = 7
THRESHOLD = 12


def run():
    s1, s2 = load_protein_pair()
    matrix = load_matrix()
    comparison = build_comparison_matrix(s1, s2, matrix)

    binary = comparison_to_dotplot(comparison, WINDOW, THRESHOLD)
    scored = comparison_to_scored_dotplot(comparison, WINDOW, THRESHOLD)

    kept = int(binary.sum())
    zeros = kept - int(np.count_nonzero(scored))
    survivors = comparison[binary > 0]

    print(f'w = {WINDOW}, threshold = {THRESHOLD}')
    print(f'  1s-and-0s dot-plot: {kept} cells on ({density(binary):.2%})')
    print(f'  dot-plot with values: the same {kept} cells, with scores '
          f'between {scored.min():.0f} and {scored.max():.0f} '
          f'({zeros} of them are worth exactly 0 and end up indistinguishable '
          f'from the background)')
    print(f'  average score of what survives the filter: '
          f'{np.mean(survivors):.2f} (against {comparison.mean():.2f} in the '
          f'full matrix)')
    print('\nThe difference from 3.2b is that now a cell that survives the '
          'filter isn\'t worth the same as any other: a tryptophan against a '
          'tryptophan (11 in BLOSUM62) weighs more than a barely conservative '
          'substitution (1). That is what the path in 3.2e takes advantage of.')

    plot_comparison_matrix(scored, s1, s2,
                           title=f'Dot-plot with real values '
                                 f'(w = {WINDOW}, threshold = {THRESHOLD})')
    return scored
