"""3.2c) Combine ahora este modulo con los realizados para ADN y realice
alineamientos de a pares entre distintas proteinas. Explore como varia el
resultado en funcion de la identidad de secuencia, el uso de distintos tipos de
filtros, y umbrales. Compare sus resultados con los obtenidos por un programa
como BLAST.
"""

from alignment import align_from_path, alignment_stats, print_alignment
from dotplot import build_comparison_matrix, comparison_to_dotplot, density
from paths import find_path
from sequences import (DISTANT_PROTEIN, data_path, load_fasta,
                       load_protein_pair)
from substitution import load_matrix

WINDOW = 7
THRESHOLD = 12
THRESHOLDS = [5, 12, 20]


def align(s1, s2, matrix, window=WINDOW, threshold=THRESHOLD):
    comparison = build_comparison_matrix(s1, s2, matrix)
    dotplot = comparison_to_dotplot(comparison, window, threshold)
    path = find_path(dotplot)
    row1, row2 = align_from_path(s1, s2, path)
    return dotplot, row1, row2


def run():
    matrix = load_matrix()
    beta, delta = load_protein_pair()
    myoglobin = load_fasta(data_path(DISTANT_PROTEIN))

    pairs = {
        'beta-globina vs. delta-globina (parecidas)': (beta, delta),
        'beta-globina vs. mioglobina (lejanas)': (beta, myoglobin),
    }

    for label, (s1, s2) in pairs.items():
        print(f'\n--- {label} ---')
        _, row1, row2 = align(s1, s2, matrix)
        print_alignment(row1, row2)

    print(f'\n--- efecto del umbral (w = {WINDOW}) ---')
    print(f'{"par":<44} {"umbral":>7} {"dots":>7} {"identidad":>10} {"gaps":>6}')
    for label, (s1, s2) in pairs.items():
        for threshold in THRESHOLDS:
            dotplot, row1, row2 = align(s1, s2, matrix, threshold=threshold)
            stats = alignment_stats(row1, row2)
            print(f'{label[:44]:<44} {threshold:>7} {density(dotplot):>7.1%} '
                  f'{stats["identity"]:>10.1%} {stats["gaps"]:>6}')

    print("""
El par parecido aguanta cualquier umbral: la diagonal esta tan poblada que el
camino la sigue igual y el alineamiento no cambia. El par lejano es el que se
rompe, pero no como uno esperaria: con umbral alto el dot-plot queda casi
vacio, el camino se queda sin casillas en que apoyarse y sale derecho por la
diagonal sin engancharse con nada. Hace MENOS gaps, no mas, y la identidad se
desploma. La senal la borro el filtro, no el algoritmo de busqueda.

Contra BLAST (blastp, "Align two or more sequences"): para beta vs. delta la
identidad deberia dar practicamente igual, porque el alineamiento es casi todo
diagonal. Para beta vs. mioglobina BLAST reporta un alineamiento local mas
corto y con mejor identidad que el nuestro, que es global; ahi la diferencia no
es un error del camino sino del modelo (global vs. local).
""")
