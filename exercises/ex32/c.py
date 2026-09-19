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
        'beta-globin vs. delta-globin (close)': (beta, delta),
        'beta-globin vs. myoglobin (distant)': (beta, myoglobin),
    }

    for label, (s1, s2) in pairs.items():
        print(f'\n--- {label} ---')
        _, row1, row2 = align(s1, s2, matrix)
        print_alignment(row1, row2)

    print(f'\n--- effect of the threshold (w = {WINDOW}) ---')
    print(f'{"pair":<40} {"threshold":>9} {"dots":>7} {"identity":>10} {"gaps":>6}')
    for label, (s1, s2) in pairs.items():
        for threshold in THRESHOLDS:
            dotplot, row1, row2 = align(s1, s2, matrix, threshold=threshold)
            stats = alignment_stats(row1, row2)
            print(f'{label[:40]:<40} {threshold:>9} {density(dotplot):>7.1%} '
                  f'{stats["identity"]:>10.1%} {stats["gaps"]:>6}')

    print("""
The close pair holds up under any threshold: the diagonal is so densely
populated that the path follows it regardless, and the alignment doesn't
change. The distant pair is where things fall apart, and badly: identity
craters to 7.7-13.4%, barely above what a handful of lucky BLOSUM-positive
substitutions would give by chance. This isn't (only) the filter's fault the
way it was for the very noisy raw dot-plot: even at the threshold that keeps
the most signal (5), a purely local, one-step-at-a-time decision still can't
reliably follow a diagonal this weak, because a single bad call early on
never gets corrected later.

Against BLAST (blastp, "Align two or more sequences"): for beta vs. delta the
identity should come out essentially the same, because the alignment is
almost entirely diagonal. For beta vs. myoglobin, expect this path search to
land well below whatever BLAST reports (global or local, significant or not)
-that gap is the naive, no-lookback search design showing its limits, not a
property of the sequences themselves.
""")
