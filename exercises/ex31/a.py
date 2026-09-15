"""3.1a) Cree una matriz de ceros (de dimension de acuerdo a las longitudes de
seq1 y seq2) y luego recorra todas las casillas poniendo 1 donde hay un match y
0 donde hay mismatch.

Consejo del TP: pruebe con una matriz mas pequena y luego vaya por la completa.
"""

from dotplot import build_dotplot, density, preview
from sequences import load_pair

SMALL = 20   # cuantas bases usa la prueba chica


def run():
    s1, s2 = load_pair()

    # Primero la matriz chica, para poder mirarla entera.
    print(f'--- prueba con las primeras {SMALL} bases ---')
    print(f's1: {s1[:SMALL]}')
    print(f's2: {s2[:SMALL]}')
    small = build_dotplot(s1[:SMALL], s2[:SMALL])
    print(f'matriz {small.shape[0]} x {small.shape[1]}, '
          f'{int(small.sum())} matches ({density(small):.1%} de las casillas)')
    print(preview(small, rows=SMALL, cols=SMALL))

    # Y ahora las secuencias completas.
    print(f'\n--- secuencias completas ---')
    print(f's1: {len(s1)} bases   s2: {len(s2)} bases')
    dotplot = build_dotplot(s1, s2)
    print(f'matriz {dotplot.shape[0]} x {dotplot.shape[1]}, '
          f'{int(dotplot.sum())} matches ({density(dotplot):.1%} de las casillas)')
    print('esquina superior izquierda:')
    print(preview(dotplot))

    print('\nCon 4 bases equiprobables, dos posiciones cualquiera coinciden 1 de '
          'cada 4 veces: el ~25% de fondo es ruido, no senal. Eso es lo que hay '
          'que filtrar en 3.1c.')
    return dotplot
