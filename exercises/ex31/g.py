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

    print('resultado de nuestro algoritmo:')
    print(f'  largo del alineamiento: {stats["length"]}')
    print(f'  identidad: {stats["identity"]:.1%} '
          f'({stats["matches"]}/{stats["length"]})')
    print(f'  gaps: {stats["gaps"]} (pasos V={counts["V"]}, H={counts["H"]})')

    paths = []
    for filename, sequence in zip(EXPORTS, (s1, s2)):
        paths.append(write_fasta(data_path(filename), sequence, filename.split('.')[0]))
    print('\nsecuencias exportadas para pegar en BLAST:')
    for exported in paths:
        print(f'  {exported}')

    print("""
Como comparar:
  1. Entrar a https://blast.ncbi.nlm.nih.gov -> Nucleotide BLAST.
  2. Tildar "Align two or more sequences" y pegar los dos FASTA exportados.
  3. En Program selection elegir "blastn" (mas sensible que megablast, que
     esta pensado para secuencias casi identicas).
  4. Comparar contra los numeros de arriba: identidad, largo del alineamiento
     y cantidad de gaps.

Que esperar: BLAST es local (alinea solo los tramos que valen la pena) y
nuestro camino es global (recorre las secuencias enteras de punta a punta), asi
que BLAST suele reportar un alineamiento mas corto y con identidad mas alta.
Si dentro del tramo que BLAST reporta nuestro camino sigue la misma diagonal,
el algoritmo de busqueda de camino esta bien.
""")
    return stats
