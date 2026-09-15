"""Del camino sobre el dot-plot al alineamiento en pantalla (3.1 h).

Formato pedido por el TP, 3 filas:

    ACGTTCAGTAG     f1: secuencia 1 y sus gaps
    AC---CACTTG     f2: secuencia 2 y sus gaps
    **    *** *     f3: indicador de match

Los gaps salen de los pasos del camino. Con la convencion de este codigo (s1 en
las filas, s2 en las columnas):

    D  avanza en las dos secuencias  -> se alinean s1[i] y s2[j]
    V  avanza solo en s1 (la vertical) -> gap en s2
    H  avanza solo en s2 (la horizontal) -> gap en s1
"""

from paths import path_steps

GAP = '-'


def align_from_path(s1, s2, path):
    """Devuelve (fila1, fila2): s1 y s2 con los gaps que impone el camino."""
    if not path:
        return '', ''

    i, j = path[0]
    # La primera casilla del camino ya es un par alineado.
    row1, row2 = s1[i], s2[j]

    for step in path_steps(path):
        if step == 'D':
            i, j = i + 1, j + 1
            row1 += s1[i]
            row2 += s2[j]
        elif step == 'V':
            i += 1
            row1 += s1[i]
            row2 += GAP
        else:                      # H
            j += 1
            row1 += GAP
            row2 += s2[j]

    return row1, row2


def match_row(row1, row2, symbol='*'):
    """Tercera fila del alineamiento: `symbol` donde las dos filas coinciden."""
    marks = ''
    for a, b in zip(row1, row2):
        marks += symbol if a == b and a != GAP else ' '
    return marks


def format_alignment(row1, row2, match_first=False, width=60):
    """Arma el bloque de 3 filas listo para imprimir, cortado en `width`.

    match_first=True intercambia f2 y f3, la otra variante que menciona el TP.
    """
    marks = match_row(row1, row2)
    blocks = []
    for start in range(0, len(row1), width):
        chunk1 = row1[start:start + width]
        chunk2 = row2[start:start + width]
        chunk_marks = marks[start:start + width]
        rows = [chunk1, chunk_marks, chunk2] if match_first else [chunk1, chunk2, chunk_marks]
        header = f'{start + 1}-{min(start + width, len(row1))}'
        blocks.append(f'{header}\n' + '\n'.join(rows))
    return '\n\n'.join(blocks)


def alignment_stats(row1, row2):
    """Resumen del alineamiento: largo, matches, mismatches, gaps e identidad."""
    matches = mismatches = gaps = 0
    for a, b in zip(row1, row2):
        if a == GAP or b == GAP:
            gaps += 1
        elif a == b:
            matches += 1
        else:
            mismatches += 1
    length = len(row1)
    return {
        'length': length,
        'matches': matches,
        'mismatches': mismatches,
        'gaps': gaps,
        'identity': matches / length if length else 0.0,
    }


def print_alignment(row1, row2, match_first=False, width=60):
    """Imprime el alineamiento y su resumen, que es lo que consumen 3.1h y 3.2c."""
    print(format_alignment(row1, row2, match_first=match_first, width=width))
    stats = alignment_stats(row1, row2)
    print(f"\nlargo del alineamiento: {stats['length']}"
          f"  |  matches: {stats['matches']}"
          f"  |  mismatches: {stats['mismatches']}"
          f"  |  gaps: {stats['gaps']}"
          f"  |  identidad: {stats['identity']:.1%}")
    return stats
