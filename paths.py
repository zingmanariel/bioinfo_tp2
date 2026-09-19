"""Path search over the dot-plot: a greedy walk (exercises 3.1 d/e, 3.2 c/e).

A path is a list of coordinates [(i, j), ...] starting at the top-left corner
and advancing to the bottom-right corner. Between two consecutive coordinates
there are only three possible steps:

    D (diagonal)   i+1, j+1   aligns s1[i] with s2[j] (match or mismatch)
    V (vertical)   i+1, j     gap in s2
    H (horizontal) i,   j+1   gap in s1

How each step is decided (the TP's "divide and conquer" hint, read as a local
decision taken one step at a time): standing at (i, j), look at the value of
the three neighboring cells -diagonal, vertical, horizontal- and move to
whichever is highest right now. That's it: there is no table of accumulated
scores for the whole matrix, and a decision already made is never revisited.
This is the opposite of dynamic programming (see dyn_align.py, which does
keep the best accumulated score at every cell): a greedy path can make a bad
call early on and stay off track for the rest of the sequence, something a DP
never does because it considers every path before committing to one.
"""

import numpy as np


def find_path(dotplot, start=(0, 0)):
    """Greedy path over the dot-plot, as a list of (i, j) coordinates.

    At each step it compares the value of the three neighboring cells
    (diagonal, vertical, horizontal) and moves to whichever is highest right
    now; ties go to the diagonal, then to vertical. Works the same on a 0/1
    dot-plot and on one with real values (3.2d): with real values the
    comparison is no longer "is there a match", it's "how good is it", so a
    conservative substitution can outweigh a weak match.
    """
    array = np.asarray(dotplot, dtype=float)
    rows, cols = array.shape
    i, j = start
    path = [(i, j)]

    while i < rows - 1 or j < cols - 1:
        options = []
        if i + 1 < rows and j + 1 < cols:
            options.append((array[i + 1][j + 1], (i + 1, j + 1)))   # D
        if i + 1 < rows:
            options.append((array[i + 1][j], (i + 1, j)))           # V
        if j + 1 < cols:
            options.append((array[i][j + 1], (i, j + 1)))           # H
        _, (i, j) = max(options, key=lambda option: option[0])
        path.append((i, j))

    return path


# ---------------------------------------------------------------------------
# Utilities over an already computed path
# ---------------------------------------------------------------------------

def path_steps(path):
    """Translate the path into its D/V/H step string, to inspect it.

    This is what alignment.align_from_path() consumes to know where the gaps go.
    """
    steps = ''
    for (i0, j0), (i1, j1) in zip(path, path[1:]):
        di, dj = i1 - i0, j1 - j0
        if (di, dj) == (1, 1):
            steps += 'D'
        elif (di, dj) == (1, 0):
            steps += 'V'
        elif (di, dj) == (0, 1):
            steps += 'H'
        else:
            raise ValueError(f'Invalid step in path: {(i0, j0)} -> {(i1, j1)}')
    return steps


def step_counts(path):
    """How many steps of each type the path has: {'D': n, 'V': n, 'H': n}."""
    steps = path_steps(path)
    return {step: steps.count(step) for step in 'DVH'}


def path_score(dotplot, path):
    """Sum of the dot-plot values along the path: a quick way to compare two
    paths over the same matrix."""
    array = np.asarray(dotplot, dtype=float)
    return float(sum(array[i][j] for i, j in path))
