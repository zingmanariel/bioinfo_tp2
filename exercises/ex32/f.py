"""3.2f) Utilice los modulos d) y e) para realizar nuevos alineamientos y compare
los resultados obtenidos con respecto a los anteriores (3.2c) y con los
obtenidos por BLAST (u otro similar).

Punto opcional del TP.
"""

from alignment import align_from_path, alignment_stats
from dotplot import (build_comparison_matrix, comparison_to_dotplot,
                     comparison_to_scored_dotplot, density)
from paths import find_path, path_score
from sequences import (DISTANT_PROTEIN, data_path, load_fasta,
                       load_protein_pair)
from substitution import load_matrix

WINDOW = 7
THRESHOLDS = [0, 5, 8, 12, 20]   # 0 = unfiltered, kept only as a reference point


def summarize(s1, s2, comparison, path):
    """Summary of a path: alignment identity and BLOSUM score.

    The score is always measured over the full comparison matrix, not over
    the filtered dot-plot: it's the only way to compare, with the same
    yardstick, two paths computed over different matrices.
    """
    stats = alignment_stats(*align_from_path(s1, s2, path))
    stats['score'] = path_score(comparison, path)
    return stats


def compare(s1, s2, comparison, threshold):
    """The same pair aligned with the 1s-and-0s dot-plot and with the real-valued one."""
    binary = comparison_to_dotplot(comparison, WINDOW, threshold)
    scored = comparison_to_scored_dotplot(comparison, WINDOW, threshold)
    return {
        'dots': density(binary),
        '1s and 0s': summarize(s1, s2, comparison, find_path(binary)),
        'real values': summarize(s1, s2, comparison, find_path(scored)),
    }


def run():
    matrix = load_matrix()
    beta, delta = load_protein_pair()
    myoglobin = load_fasta(data_path(DISTANT_PROTEIN))

    pairs = {
        'beta vs. delta (close)': (beta, delta),
        'beta vs. myoglobin (distant)': (beta, myoglobin),
    }

    print(f'w = {WINDOW}   (threshold 0 = unfiltered comparison matrix)')

    for label, (s1, s2) in pairs.items():
        comparison = build_comparison_matrix(s1, s2, matrix)
        print(f'\n--- {label} ---')
        print(f'{"thresh":>6}{"dots":>8}   '
              f'{"1s and 0s (3.2c)":>28}   {"real values (3.2e)":>28}')
        print(f'{"":>6}{"":>8}   {"id":>8}{"gaps":>6}{"score":>9}   '
              f'{"":>5}{"id":>8}{"gaps":>6}{"score":>9}')
        for threshold in THRESHOLDS:
            result = compare(s1, s2, comparison, threshold)
            binary, scored = result['1s and 0s'], result['real values']
            print(f'{threshold:>6}{result["dots"]:>8.2%}   '
                  f'{binary["identity"]:>8.1%}{binary["gaps"]:>6}'
                  f'{binary["score"]:>9.0f}   '
                  f'{"":>5}{scored["identity"]:>8.1%}{scored["gaps"]:>6}'
                  f'{scored["score"]:>9.0f}')

    print("""
Reading the table:

1. The close pair is not perfectly stable either. The binary path stays at
   93.2%/0 gaps for every threshold, but the real-valued path drifts (78-90%
   identity, 6-8 gaps depending on the threshold) and at thresholds 8+ it
   actually scores higher in raw BLOSUM (736 vs. 727) while showing lower
   identity: it's optimizing total substitution score, not match count, so it
   will trade a couple of matches for a pair of gaps if that raises the
   score.

2. The distant pair does not show the real-valued path winning. Unfiltered
   (threshold 0) it's clearly worse than the binary one: 51 gaps and a very
   negative score (-45) against 23 gaps and +16 for the binary path. Handing
   the greedy rule real magnitudes without any filtering just gives it more
   ways to talk itself into a gap chasing a marginally-less-bad neighbor, and
   it wanders more, not less. From threshold 8 up, both variants converge to
   the same low floor (7.7% at 12 and 20): once the filter has stripped out
   almost everything, there's barely any surviving structure left for either
   version to disagree about.

3. So for this algorithm, filtering and using real values are not a package
   deal that reinforces itself the way one might expect: on the close pair,
   real values only add noise to an already-solved case; on the distant pair,
   real values without filtering make things worse, and with heavy filtering
   both variants end up in the same bad place anyway.

Careful reading the identity column as if it were the alignment's grade: the
real-valued path is maximizing BLOSUM score, not number of matches. It can
prefer a conservative substitution (I for V, K for R) over an isolated match
surrounded by expensive substitutions, and that is exactly what a real
scoring-based aligner does -it's just that a purely local, one-step decision
rule doesn't use that extra information reliably.

Against BLAST: for the close pair, expect all of these alignments (3.2c, 3.2e
and BLAST) to land close to each other, since the alignment is almost
entirely diagonal either way. For the distant pair, expect this path search
to land well below BLAST regardless of filter or threshold -that gap comes
from the search itself being local and one-step-at-a-time, not from any
detail of how the dot-plot was filtered.
""")
