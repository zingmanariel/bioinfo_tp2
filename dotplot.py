"""Dot-plots: la matriz de comparacion entre dos secuencias y los filtros que
la hacen legible.

Ejercicio 3.1 a) y c) para ADN, y 3.2 a), b) y d) para proteinas.
Todo lo que se programa en este modulo es una matriz: numpy hace el resto.

Convencion que usa todo el codigo: s1 va en las FILAS (eje vertical, indice i)
y s2 en las COLUMNAS (eje horizontal, indice j).
"""

import numpy as np


# ---------------------------------------------------------------------------
# 3.1a) Dot-plot de ADN: matriz de 1s (match) y 0s (mismatch)
# ---------------------------------------------------------------------------

def build_dotplot(s1, s2):
    """Matriz len(s1) x len(s2) con 1 donde s1[i] == s2[j] y 0 si no."""
    dotplot = np.zeros((len(s1), len(s2)))
    for i, a in enumerate(s1):
        for j, b in enumerate(s2):
            if a == b:
                dotplot[i][j] = 1
    return dotplot


# ---------------------------------------------------------------------------
# 3.1c) Filtros por ventana
# ---------------------------------------------------------------------------
#
# El grafico del TP muestra decisiones independientes, que aca son parametros
# de la misma funcion:
#
#   window       tamano de la ventana (w x w)
#   threshold    umbral: cuantos matches tiene que haber dentro de la ventana
#                para que "prenda" (umbral alto vs. umbral bajo)
#   overlapping  ventanas solapadas (paso 1) o adyacentes (paso w)
#   reduce       el dot-plot se REEMPLAZA por el de las ventanas (matriz mas
#                chica, "filtro y reduccion") o las ventanas se usan solo como
#                filtro y la matriz mantiene su tamano ("filtro sin reduccion")
#   diagonal     que se cuenta adentro de la ventana: todos los matches del
#                bloque w x w (lo que dibuja el TP) o solo los de su diagonal
#                (lo que hace un dot-plot clasico; ver la nota de abajo)
#
# Nota sobre `diagonal`: contando todo el bloque, una ventana sobre la diagonal
# verdadera tiene w matches reales + ruido, y una ventana cualquiera de ADN
# tiene en promedio w*w/4 matches de fondo. Para w=10 eso es 10 + 22 contra 25:
# casi no hay contraste, y por eso el filtro de bloque limpia poco. Contando
# solo la diagonal de la ventana se comparan w matches contra w/4, que es la
# unica forma de que el filtro realmente descongestione. El ejercicio 3.1c
# corre las dos variantes y muestra los numeros.
# ---------------------------------------------------------------------------

def _window_starts(size, window, overlapping):
    """Coordenada donde arranca cada ventana sobre un eje de largo `size`."""
    if overlapping:
        return list(range(0, max(size - window + 1, 1)))
    return list(range(0, size, window))


def _window_count(block, diagonal):
    """Cuantos matches tiene la ventana: los de su diagonal o los del bloque."""
    if diagonal:
        return float(np.trace(block))
    return float(block.sum())


def filter_dotplot(dotplot, window, threshold, overlapping=True, reduce=False,
                   diagonal=False):
    """Aplica el filtro de ventana sobre un dot-plot de 0s y 1s.

    Devuelve otra matriz: del mismo tamano si reduce=False (se conservan los
    dots originales de las ventanas que superan el umbral y se apaga el resto),
    o una matriz de una celda por ventana si reduce=True.
    """
    array = np.asarray(dotplot, dtype=float)
    rows, cols = array.shape
    row_starts = _window_starts(rows, window, overlapping)
    col_starts = _window_starts(cols, window, overlapping)

    if reduce:
        result = np.zeros((len(row_starts), len(col_starts)))
    else:
        result = np.zeros_like(array)

    for r, i in enumerate(row_starts):
        for c, j in enumerate(col_starts):
            block = array[i:i + window, j:j + window]
            if _window_count(block, diagonal) < threshold:
                continue
            if reduce:
                result[r][c] = 1
            elif diagonal:
                # solo prende la diagonal de la ventana, que es lo que se conto
                offsets = np.arange(min(block.shape))
                result[i + offsets, j + offsets] = array[i + offsets, j + offsets]
            else:
                result[i:i + window, j:j + window] = block
    return result


def filter_variants(dotplot, window, threshold, diagonal=False):
    """Las cuatro combinaciones del grafico del TP, para compararlas de una.

    Devuelve un dict {etiqueta: matriz filtrada}.
    """
    variants = {}
    for overlapping in (True, False):
        for reduce in (False, True):
            label = 'solapadas' if overlapping else 'adyacentes'
            label += ', con reduccion' if reduce else ', sin reduccion'
            variants[label] = filter_dotplot(dotplot, window, threshold,
                                             overlapping=overlapping,
                                             reduce=reduce,
                                             diagonal=diagonal)
    return variants


# ---------------------------------------------------------------------------
# 3.2a) Matriz de comparacion de proteinas (puntajes de la matriz de sustitucion)
# ---------------------------------------------------------------------------

def build_comparison_matrix(s1, s2, matrix):
    """Matriz len(s1) x len(s2) con el puntaje de sustitucion de cada par.

    Igual que build_dotplot(), pero en vez de 0/1 guarda matrix[a][b]: numeros,
    positivos para sustituciones conservativas y negativos para las otras.
    `matrix` es la que devuelve substitution.load_matrix().
    """
    from substitution import score   # import local: evita un ciclo de imports

    comparison = np.zeros((len(s1), len(s2)))
    for i, a in enumerate(s1):
        for j, b in enumerate(s2):
            comparison[i][j] = score(matrix, a, b)
    return comparison


# ---------------------------------------------------------------------------
# 3.2b) y 3.2d) De la matriz de comparacion al dot-plot
# ---------------------------------------------------------------------------

def _window_mask(comparison, window, threshold, overlapping=True, diagonal=True):
    """Matriz de booleanos: True donde la ventana de la celda supera el umbral.

    overlapping=True: una ventana centrada en cada celda (paso 1), asi el
    dot-plot queda alineado con las secuencias. overlapping=False: ventanas
    adyacentes, donde prende o se apaga el bloque entero.

    diagonal=True suma los puntajes a lo largo de la diagonal de la ventana
    (una corrida de residuos parecidos), que es lo que hace legible el dot-plot
    de proteinas; diagonal=False suma todo el bloque.
    """
    array = np.asarray(comparison, dtype=float)
    rows, cols = array.shape
    mask = np.zeros((rows, cols), dtype=bool)
    half = window // 2

    if overlapping:
        for i in range(rows):
            for j in range(cols):
                if diagonal:
                    total = 0.0
                    for k in range(-half, window - half):
                        if 0 <= i + k < rows and 0 <= j + k < cols:
                            total += array[i + k][j + k]
                else:
                    block = array[max(i - half, 0):i + window - half,
                                  max(j - half, 0):j + window - half]
                    total = float(block.sum())
                mask[i][j] = total >= threshold
        return mask

    for i in _window_starts(rows, window, overlapping=False):
        for j in _window_starts(cols, window, overlapping=False):
            block = array[i:i + window, j:j + window]
            total = float(np.trace(block)) if diagonal else float(block.sum())
            if total >= threshold:
                mask[i:i + window, j:j + window] = True
    return mask


def comparison_to_dotplot(comparison, window, threshold, overlapping=True,
                          diagonal=True):
    """Convierte la matriz de puntajes en un dot-plot de 1s y 0s.

    Parametros que pide el TP: un window-size y un umbral de puntaje. La idea
    es la misma que en filter_dotplot(), pero lo que se acumula dentro de la
    ventana son puntajes y no cantidad de matches.
    """
    return _window_mask(comparison, window, threshold, overlapping,
                        diagonal).astype(float)


def comparison_to_scored_dotplot(comparison, window, threshold, overlapping=True,
                                 diagonal=True):
    """Igual que comparison_to_dotplot() pero conservando los valores REALES.

    Las celdas cuya ventana no supera el umbral quedan en 0; las demas
    mantienen su puntaje. Es el dot-plot con el que trabaja
    paths.find_path_scored() (3.2e).
    """
    mask = _window_mask(comparison, window, threshold, overlapping, diagonal)
    return np.asarray(comparison, dtype=float) * mask


# ---------------------------------------------------------------------------
# Utilidades
# ---------------------------------------------------------------------------

def density(dotplot):
    """Fraccion de celdas prendidas: sirve para ver cuanto "descongestiona" un
    filtro antes de mirar el grafico."""
    array = np.asarray(dotplot)
    return float(np.count_nonzero(array)) / array.size


def preview(dotplot, rows=20, cols=60, on='#', off='.'):
    """Esquina superior izquierda del dot-plot como texto, para mirarlo sin
    graficar (el "pruebe con una matriz mas pequena" del TP)."""
    array = np.asarray(dotplot)
    lines = []
    for row in array[:rows, :cols]:
        lines.append(''.join(on if value else off for value in row))
    return '\n'.join(lines)
