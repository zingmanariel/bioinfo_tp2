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

    print(f'--- primeros {BLOCKS_SHOWN} bloques de nuestro alineamiento ---')
    blocks = format_alignment(row1, row2, width=WIDTH).split('\n\n')
    print('\n\n'.join(blocks[:BLOCKS_SHOWN]))

    print(f'\nnuestro alineamiento: largo {stats["length"]}, '
          f'identidad {stats["identity"]:.1%}, gaps {stats["gaps"]}')

    for filename, sequence in zip(EXPORTS, (s1, s2)):
        write_fasta(data_path(filename), sequence, filename.split('.')[0])
    print(f'secuencias exportadas a data/{EXPORTS[0]} y data/{EXPORTS[1]}')

    print("""
Que mirar al comparar con BLAST (bl2seq, blastn):
  - Extension: BLAST alinea por tramos (HSPs) locales; nuestro camino es
    global, asi que cubre las secuencias enteras, incluidas las puntas donde no
    hay homologia (UTRs). Esas puntas son las que bajan nuestra identidad.
  - Identidad: dentro del tramo que BLAST reporta deberia coincidir con la
    nuestra; si difiere mucho, el camino se esta yendo de la diagonal.
  - Gaps: BLAST penaliza abrir un gap mas caro que extenderlo (affine gaps);
    nuestro modelo cobra lo mismo por cada paso V/H (gap lineal), asi que
    tendemos a repartir mas gaps cortos donde BLAST pone uno solo largo.
""")
    return stats
