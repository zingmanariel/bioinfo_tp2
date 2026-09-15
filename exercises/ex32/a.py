"""3.2a) Busque una matriz de sustitucion (Blossum60?). Levantela en python y
luego de analizarla programe un codigo que haga un dot-plot para el
alineamiento de 2 secuencias proteicas con la misma.

El problema que aparece: la matriz de comparacion esta llena de numeros, y lo
que hace falta para el camino es una matriz de 1s y 0s. Esa conversion es 3.2b.
"""

import numpy as np

from dotplot import build_comparison_matrix
from plots import plot_comparison_matrix
from sequences import load_protein_pair
from substitution import DEFAULT_MATRIX, available, describe, load_matrix


def run():
    print(f'matrices disponibles en Biopython:\n  {", ".join(available())}\n')
    matrix = load_matrix(DEFAULT_MATRIX)
    describe(matrix, DEFAULT_MATRIX)

    s1, s2 = load_protein_pair()
    print(f'\nproteinas: {len(s1)} y {len(s2)} residuos')
    comparison = build_comparison_matrix(s1, s2, matrix)
    print(f'matriz de comparacion {comparison.shape[0]} x {comparison.shape[1]}: '
          f'puntajes entre {comparison.min():.0f} y {comparison.max():.0f}, '
          f'promedio {comparison.mean():.2f}')
    positive = float(np.count_nonzero(comparison > 0)) / comparison.size
    print(f'{positive:.1%} de las casillas tiene puntaje positivo; ese es el '
          f'"dot-plot ingenuo" (prender donde el puntaje es > 0) y sale tan '
          f'congestionado como el de ADN sin filtrar.')

    plot_comparison_matrix(comparison, s1, s2,
                           title=f'Matriz de comparacion ({DEFAULT_MATRIX})')
    return comparison
