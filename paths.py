"""Busqueda del camino sobre el dot-plot.

Ejercicios 3.1 d) y e) para el dot-plot de 0s y 1s, y 3.2 e) para el dot-plot
con valores reales.

Un camino es una lista de coordenadas [(i, j), ...] que arranca arriba a la
izquierda y avanza hacia abajo a la derecha. Entre dos coordenadas consecutivas
solo hay tres pasos posibles:

    D (diagonal)   i+1, j+1   alinea s1[i] con s2[j] (match o mismatch)
    V (vertical)   i+1, j     gap en s2
    H (horizontal) i,   j+1   gap en s1

Como se decide cada paso (el consejo de "dividir y conquistar" del TP): el
mejor camino que llega a la casilla (i, j) es el mejor de los tres caminos que
llegan a sus vecinas (i-1, j-1), (i-1, j) e (i, j-1), mas lo que suma el paso
final. Es decir, el problema grande se arma con la solucion de los problemas
chicos, y alcanza con recorrer la matriz una sola vez guardando, en cada
casilla, el mejor puntaje y de donde vino. Al final se vuelve hacia atras
desde la mejor casilla del borde derecho o inferior: eso es el camino.
"""

import numpy as np


# ---------------------------------------------------------------------------
# 3.1d) Recorrer el dot-plot de 0s y 1s
# ---------------------------------------------------------------------------

def find_path(dotplot, start=(0, 0), gap_penalty=-1.0):
    """Mejor camino sobre el dot-plot, como lista de coordenadas (i, j).

    Cada paso D suma el valor de la casilla (1 si es match, 0 si no) y cada
    paso V/H cuesta gap_penalty, asi que el camino prefiere la diagonal
    mientras haya matches y solo mete un gap cuando le conviene.
    """
    return _best_path(dotplot, start, gap_penalty)


# ---------------------------------------------------------------------------
# 3.2e) El mismo recorrido, pero sobre un dot-plot con valores reales
# ---------------------------------------------------------------------------

def find_path_scored(scored_dotplot, start=(0, 0), gap_penalty=-4.0):
    """Camino sobre el dot-plot de valores reales de 3.2d).

    Es el mismo algoritmo que find_path(): lo unico que cambia es que ahora
    cada casilla vale su puntaje de sustitucion en vez de 0 o 1, y que el gap
    cuesta lo que cuesta en esa escala (-4 es el orden de un gap en BLOSUM62).
    """
    return _best_path(scored_dotplot, start, gap_penalty)


def _best_path(dotplot, start, gap_penalty):
    """Programacion dinamica sobre la matriz + reconstruccion del camino."""
    i0, j0 = start
    values = np.asarray(dotplot, dtype=float)[i0:, j0:].tolist()
    if not values or not values[0]:
        return []
    rows, cols = len(values), len(values[0])

    # score[i][j]: puntaje del mejor camino que va de (0, 0) a (i, j)
    # move[i][j]:  con que paso se llego a (i, j) ('S' = arranque)
    score = [[0.0] * cols for _ in range(rows)]
    move = [[''] * cols for _ in range(rows)]
    score[0][0] = values[0][0]
    move[0][0] = 'S'

    for i in range(rows):
        for j in range(cols):
            if i == 0 and j == 0:
                continue
            options = []
            if i > 0 and j > 0:
                options.append((score[i - 1][j - 1] + values[i][j], 'D'))
            if i > 0:
                options.append((score[i - 1][j] + gap_penalty, 'V'))
            if j > 0:
                options.append((score[i][j - 1] + gap_penalty, 'H'))
            # max() se queda con el primero de los empatados, y D va primero:
            # ante igualdad de puntaje conviene el paso diagonal.
            score[i][j], move[i][j] = max(options, key=lambda option: option[0])

    return _traceback(score, move, i0, j0)


def _best_end(score):
    """Mejor casilla del borde derecho o inferior: ahi termina el camino.

    El TP pide llegar al limite derecho o inferior (idealmente a la esquina),
    asi que el final no se fuerza en (rows-1, cols-1): se elige el borde que
    mejor puntaje da.
    """
    rows, cols = len(score), len(score[0])
    best = (rows - 1, cols - 1)
    for i in range(rows):
        if score[i][cols - 1] > score[best[0]][best[1]]:
            best = (i, cols - 1)
    for j in range(cols):
        if score[rows - 1][j] > score[best[0]][best[1]]:
            best = (rows - 1, j)
    return best


def _traceback(score, move, i0, j0):
    """Rehace el camino hacia atras desde el final hasta el arranque."""
    i, j = _best_end(score)
    path = []
    while True:
        path.append((i + i0, j + j0))
        step = move[i][j]
        if step == 'S':
            break
        if step == 'D':
            i, j = i - 1, j - 1
        elif step == 'V':
            i -= 1
        else:
            j -= 1
    path.reverse()
    return path


# ---------------------------------------------------------------------------
# Utilidades sobre un camino ya calculado
# ---------------------------------------------------------------------------

def path_steps(path):
    """Traduce el camino a la cadena de pasos D/V/H, para inspeccionarlo.

    Es lo que consume alignment.align_from_path() para saber donde van los gaps.
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
            raise ValueError(f'Paso invalido en el camino: {(i0, j0)} -> {(i1, j1)}')
    return steps


def step_counts(path):
    """Cuantos pasos de cada tipo tiene el camino: {'D': n, 'V': n, 'H': n}."""
    steps = path_steps(path)
    return {step: steps.count(step) for step in 'DVH'}


def path_score(dotplot, path):
    """Suma de los valores del dot-plot a lo largo del camino: una medida
    rapida para comparar dos caminos sobre la misma matriz."""
    array = np.asarray(dotplot, dtype=float)
    return float(sum(array[i][j] for i, j in path))
