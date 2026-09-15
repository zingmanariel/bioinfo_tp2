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

    print(f'w = {WINDOW}, umbral = {THRESHOLD}')
    print(f'  dot-plot de 1s y 0s: {kept} casillas prendidas ({density(binary):.2%})')
    print(f'  dot-plot con valores: las mismas {kept} casillas, con puntajes '
          f'entre {scored.min():.0f} y {scored.max():.0f} '
          f'({zeros} de ellas valen exactamente 0 y quedan indistinguibles '
          f'del fondo)')
    print(f'  puntaje promedio de lo que sobrevive al filtro: '
          f'{np.mean(survivors):.2f} (contra {comparison.mean():.2f} en la '
          f'matriz completa)')
    print('\nLa diferencia con 3.2b es que ahora una casilla que sobrevive al '
          'filtro no vale lo mismo que cualquier otra: un triptofano contra un '
          'triptofano (11 en BLOSUM62) pesa mas que una sustitucion apenas '
          'conservativa (1). Eso es lo que aprovecha el camino de 3.2e.')

    plot_comparison_matrix(scored, s1, s2,
                           title=f'Dot-plot con valores reales '
                                 f'(w = {WINDOW}, umbral = {THRESHOLD})')
    return scored
