"""3.1f) Analice sobre el camino el efecto de utilizar secuencias con diferente
grado de identidad, diferente largo, y el uso de los filtros programados en el
punto anterior.

Para poder mover la identidad y el largo a voluntad, los pares de este
ejercicio se generan con sequences.mutate() y sequences.insert_indels() en vez
de leerse de data/: asi se sabe de antemano cuanta identidad tiene el par y se
puede ver que tan bien la recupera el camino.
"""

from alignment import align_from_path, alignment_stats
from dotplot import build_dotplot, density, filter_dotplot
from paths import find_path, step_counts
from plots import plot_dotplot_grid
from sequences import (generate_random_nt_sequence, insert_indels, mutate)

LENGTH = 300
IDENTITIES = [0.95, 0.85, 0.70, 0.50]
INDELS = 6
WINDOW = 10
THRESHOLD = 6


def measure(s1, s2, dotplot):
    """Run the path search over the dot-plot and return the alignment summary."""
    path = find_path(dotplot)
    row1, row2 = align_from_path(s1, s2, path)
    stats = alignment_stats(row1, row2)
    stats['gap_steps'] = step_counts(path)['V'] + step_counts(path)['H']
    return path, stats


def row(label, dotplot, stats):
    print(f'{label:<22} {density(dotplot):>7.2%} {stats["identity"]:>11.1%} '
          f'{stats["length"]:>8} {stats["gap_steps"]:>7}')


def header(title):
    print(f'\n--- {title} ---')
    print(f'{"":<22} {"dots":>7} {"identity":>11} {"length":>8} {"gaps":>7}')


def run():
    # 1) identity: same length, increasing mutations
    header(f'effect of identity (length {LENGTH}, no filter)')
    dotplots, paths = {}, {}
    for identity in IDENTITIES:
        s1 = generate_random_nt_sequence(LENGTH)
        s2 = mutate(s1, identity)
        dotplot = build_dotplot(s1, s2)
        path, stats = measure(s1, s2, dotplot)
        row(f'identity {identity:.0%}', dotplot, stats)
        label = f'identity {identity:.0%}'
        dotplots[label], paths[label] = dotplot, path

    # 2) different length: same mutations, plus indels
    header(f'effect of a length mismatch ({INDELS} indels)')
    for identity in IDENTITIES[:2]:
        s1 = generate_random_nt_sequence(LENGTH)
        s2 = insert_indels(mutate(s1, identity), INDELS)
        dotplot = build_dotplot(s1, s2)
        path, stats = measure(s1, s2, dotplot)
        row(f'identity {identity:.0%} + indels', dotplot, stats)

    # 3) filter: the same pair, filtered and unfiltered
    header(f'effect of the filter (w = {WINDOW}, threshold = {THRESHOLD})')
    for identity in IDENTITIES:
        s1 = generate_random_nt_sequence(LENGTH)
        s2 = mutate(s1, identity)
        dotplot = build_dotplot(s1, s2)
        filtered = filter_dotplot(dotplot, WINDOW, THRESHOLD, diagonal=True)
        _, stats = measure(s1, s2, dotplot)
        row(f'identity {identity:.0%} raw', dotplot, stats)
        _, filtered_stats = measure(s1, s2, filtered)
        row(f'identity {identity:.0%} filtered', filtered, filtered_stats)

    print("""
As true identity drops, the path opens more and more gaps instead of just
running straight through the noise: 24 gaps at 95% identity, 106 at 50%. The
greedy rule looks one step ahead at a time, so a background match on a
neighboring cell (about 25% of cells are a "match" by chance alone, since DNA
only has 4 letters) is enough to pull it off the true diagonal; measured
identity ends up well below the real one (24.1% at 50% true identity, barely
above the 25% background rate) because the path is not just "measuring"
identity anymore, it's actively chasing noise.

Indels hurt this algorithm even more than substitutions. 6 indels on top of
95%-identity sequences alone crash measured identity to 36.9% with about 100
gaps: once an indel shifts the true diagonal, the one-step lookahead has no
way to relocate it except by stumbling into a nearby match, and usually it
doesn't recover cleanly for a long stretch.

The filter is not cosmetic for this algorithm, it's load-bearing: at every
identity level it cuts gaps close to zero and pulls measured identity back
near the true value (95%: 73.4% -> 96.3% identity, gaps 40 -> 0; 70%: 26.0% ->
63.9%, gaps 130 -> 10). Removing the background noise before walking is what
lets a purely local, no-lookback decision rule work at all.
""")

    plot_dotplot_grid(dotplots, title=f'Path vs. identity (length {LENGTH})',
                      paths=paths)
