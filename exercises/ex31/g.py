"""3.1g) Busque algun programa on-line (BLAST?) que haga alineamientos y
compare los resultados de su algoritmo con el obtenido con BLAST, hasta estar
conforme de que su algoritmo de busqueda de camino es lo suficientemente bueno.

Este punto se resuelve a mano: el codigo deja las secuencias exportadas y los
numeros propios impresos, listos para contrastar contra
https://blast.ncbi.nlm.nih.gov (opcion "Align two or more sequences", bl2seq).
"""

from alignment import align_from_path, alignment_stats
from dotplot import build_dotplot
from paths import find_path, step_counts
from sequences import data_path, load_pair, write_fasta

EXPORTS = ('blast_seq1.fasta', 'blast_seq2.fasta')


def run():
    s1, s2 = load_pair()
    dotplot = build_dotplot(s1, s2)
    path = find_path(dotplot)
    row1, row2 = align_from_path(s1, s2, path)
    stats = alignment_stats(row1, row2)
    counts = step_counts(path)

    print('our path search result:')
    print(f'  alignment length: {stats["length"]}')
    print(f'  identity: {stats["identity"]:.1%} '
          f'({stats["matches"]}/{stats["length"]})')
    print(f'  gaps: {stats["gaps"]} (steps V={counts["V"]}, H={counts["H"]})')

    paths = []
    for filename, sequence in zip(EXPORTS, (s1, s2)):
        paths.append(write_fasta(data_path(filename), sequence, filename.split('.')[0]))
    print('\nsequences exported to paste into BLAST:')
    for exported in paths:
        print(f'  {exported}')

    print("""
How to compare:
  1. Go to https://blast.ncbi.nlm.nih.gov -> Nucleotide BLAST.
  2. Check "Align two or more sequences" and paste the two exported FASTA.
  3. Under Program selection pick "blastn" (more sensitive than megablast,
     which is meant for near-identical sequences).
  4. Compare against the numbers above: identity, alignment length and gaps.

What to expect: the path above comes from a purely local, one-step-at-a-time
decision (see paths.py) with no accumulated score and no way to go back on a
bad call. On real sequences the background noise (about 25% of cell pairs
match by chance alone) is enough to pull it off the true diagonal over and
over, so it can land far below both BLAST's result and the alignment
dyn_align.py's dynamic programming finds on the same pair. If it lands close
to BLAST's identity within the region BLAST reports, the path search did
well; if it's far off, that's the naive, no-lookback design showing its
limits, not a bug to chase.
""")
    return stats
