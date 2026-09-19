"""3.1i) Busque algun programa on-line (BLAST?) que haga alineamientos y compare
los resultados de su algoritmo con el obtenido con BLAST.

Igual que 3.1g, pero comparando ya el alineamiento en formato de 3 filas contra
el que devuelve BLAST: identidad, gaps y extension del alineamiento.
"""

from alignment import align_from_path, format_alignment, alignment_stats
from dotplot import build_dotplot
from paths import find_path
from sequences import data_path, load_pair, write_fasta

EXPORTS = ('blast_seq1.fasta', 'blast_seq2.fasta')
BLOCKS_SHOWN = 3
WIDTH = 60


def run():
    s1, s2 = load_pair()
    path = find_path(build_dotplot(s1, s2))
    row1, row2 = align_from_path(s1, s2, path)
    stats = alignment_stats(row1, row2)

    print(f'--- first {BLOCKS_SHOWN} blocks of our alignment ---')
    blocks = format_alignment(row1, row2, width=WIDTH).split('\n\n')
    print('\n\n'.join(blocks[:BLOCKS_SHOWN]))

    print(f'\nour alignment: length {stats["length"]}, '
          f'identity {stats["identity"]:.1%}, gaps {stats["gaps"]}')

    for filename, sequence in zip(EXPORTS, (s1, s2)):
        write_fasta(data_path(filename), sequence, filename.split('.')[0])
    print(f'sequences exported to data/{EXPORTS[0]} and data/{EXPORTS[1]}')

    print("""
What to look at when comparing with BLAST (bl2seq, blastn):
  - Coverage: BLAST aligns local stretches (HSPs); our path search here is
    global, so it covers the whole sequences, including the ends where there
    is no homology (UTRs). That is what drags our identity down.
  - Identity: on top of that, the path itself is a purely local, one-step
    decision with no accumulated score (see paths.py) -it can wander off the
    true diagonal well before reaching those ends. If it lands far below
    BLAST's identity even within the region BLAST reports, that is the greedy
    design's own limitation, not (only) the global-vs-local difference.
  - Gaps: BLAST charges more for opening a gap than for extending it (affine
    gaps); our path pays the same for every V/H step (linear gap), and on top
    of that it has no global view to decide where a single long gap would pay
    off better than several short ones.
""")
    return stats
