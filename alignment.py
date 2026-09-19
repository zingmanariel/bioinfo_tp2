"""From the path over the dot-plot to the on-screen alignment (3.1 h).

Format requested by the TP, 3 rows:

    ACGTTCAGTAG     row 1: sequence 1 and its gaps
    AC---CACTTG     row 2: sequence 2 and its gaps
    **    *** *     row 3: match indicator

The gaps come from the path's steps. With this code's convention (s1 in the
rows, s2 in the columns):

    D  advances in both sequences      -> aligns s1[i] and s2[j]
    V  advances only in s1 (vertical)   -> gap in s2
    H  advances only in s2 (horizontal) -> gap in s1
"""

from paths import path_steps

GAP = '-'


def align_from_path(s1, s2, path):
    """Returns (row1, row2): s1 and s2 with the gaps the path imposes."""
    if not path:
        return '', ''

    i, j = path[0]
    # The path's first cell is already an aligned pair.
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
    """Third row of the alignment: `symbol` where the two rows agree."""
    marks = ''
    for a, b in zip(row1, row2):
        marks += symbol if a == b and a != GAP else ' '
    return marks


def format_alignment(row1, row2, match_first=False, width=60):
    """Builds the 3-row block ready to print, wrapped at `width`.

    match_first=True swaps rows 2 and 3, the other variant the TP mentions.
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
    """Alignment summary: length, matches, mismatches, gaps and identity."""
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
    """Prints the alignment and its summary, which is what 3.1h and 3.2c consume."""
    print(format_alignment(row1, row2, match_first=match_first, width=width))
    stats = alignment_stats(row1, row2)
    print(f"\nalignment length: {stats['length']}"
          f"  |  matches: {stats['matches']}"
          f"  |  mismatches: {stats['mismatches']}"
          f"  |  gaps: {stats['gaps']}"
          f"  |  identity: {stats['identity']:.1%}")
    return stats
