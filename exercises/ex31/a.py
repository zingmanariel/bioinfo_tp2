"""3.1a) Cree una matriz de ceros (de dimension de acuerdo a las longitudes de
seq1 y seq2) y luego recorra todas las casillas poniendo 1 donde hay un match y
0 donde hay mismatch.

Consejo del TP: pruebe con una matriz mas pequena y luego vaya por la completa.
"""

from dotplot import build_dotplot, density, preview
from sequences import load_pair

SMALL = 20   # how many bases the small test uses


def run():
    s1, s2 = load_pair()

    # First the small matrix, so it can be looked at in full.
    print(f'--- test with the first {SMALL} bases ---')
    print(f's1: {s1[:SMALL]}')
    print(f's2: {s2[:SMALL]}')
    small = build_dotplot(s1[:SMALL], s2[:SMALL])
    print(f'matrix {small.shape[0]} x {small.shape[1]}, '
          f'{int(small.sum())} matches ({density(small):.1%} of the cells)')
    print(preview(small, rows=SMALL, cols=SMALL))

    # And now the full sequences.
    print(f'\n--- full sequences ---')
    print(f's1: {len(s1)} bases   s2: {len(s2)} bases')
    dotplot = build_dotplot(s1, s2)
    print(f'matrix {dotplot.shape[0]} x {dotplot.shape[1]}, '
          f'{int(dotplot.sum())} matches ({density(dotplot):.1%} of the cells)')
    print('top-left corner:')
    print(preview(dotplot))

    print('\nWith 4 equally likely bases, any two positions coincide 1 in 4 '
          'times: the ~25% background is noise, not signal. That is what has '
          'to be filtered out in 3.1c.')
    return dotplot
