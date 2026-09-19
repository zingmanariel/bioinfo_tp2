"""3.1h) Escriba un codigo que tome como entrada s1, s2 y el camino del dot-plot
y sea capaz de mostrar el alineamiento en el formato deseado:

    ACGTTCAGTAG
    AC---CACTTG
    **    *** *

Nota del TP: la dificultad esta en como poner los gaps. Ir por la Horizontal
implica poner un gap en la secuencia vertical y vice-versa. Conviene resolver
la logica en papel antes de programarla (esta escrita en alignment.py).
"""

from alignment import align_from_path, print_alignment
from dotplot import build_dotplot
from paths import find_path, path_steps
from sequences import load_pair

# The example from the assignment, to check the format on something small and verifiable.
EXAMPLE = ('ACGTTCAGTAG', 'ACCACTTG')


def align(s1, s2):
    dotplot = build_dotplot(s1, s2)
    path = find_path(dotplot)
    return path, align_from_path(s1, s2, path)


def run():
    print('--- short example from the assignment ---')
    s1, s2 = EXAMPLE
    path, (row1, row2) = align(s1, s2)
    print(f'path: {path_steps(path)}\n')
    print_alignment(row1, row2)

    print('\n--- the two sequences from data/ ---')
    s1, s2 = load_pair()
    path, (row1, row2) = align(s1, s2)
    stats = print_alignment(row1, row2)

    print('\n(the variant with the match row in the middle is obtained with '
          'print_alignment(row1, row2, match_first=True))')
    return stats
