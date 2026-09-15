"""Graficos compartidos: el dot-plot, el camino dibujado encima y la matriz de
comparacion de proteinas.

Los graficos que usa un solo ejercicio viven en el modulo de ese ejercicio.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

BLUE = '#4a6fa5'
ORANGE = '#d1852f'
GREY = '#b0b0b0'
GRID = '#e5e5e5'

# 0 = mismatch (fondo), 1 = match (punto)
DOT_CMAP = ListedColormap(['white', BLUE])

# Hasta este largo se escriben las letras de las secuencias en los ejes.
MAX_TICK_LABELS = 50


def _label_axes(ax, dotplot, s1, s2):
    """Ejes: las letras si las secuencias son cortas, los indices si no."""
    rows, cols = np.asarray(dotplot).shape

    if s1 is not None and len(s1) == rows and rows <= MAX_TICK_LABELS:
        ax.set_yticks(range(rows))
        ax.set_yticklabels(list(s1), fontsize=7)
    if s2 is not None and len(s2) == cols and cols <= MAX_TICK_LABELS:
        ax.set_xticks(range(cols))
        ax.set_xticklabels(list(s2), fontsize=7)

    ax.set_ylabel('s1' if s1 is None else f's1 ({rows})')
    ax.set_xlabel('s2' if s2 is None else f's2 ({cols})')


def _draw_dotplot(ax, dotplot, s1=None, s2=None, title=None):
    array = np.asarray(dotplot, dtype=float)
    ax.imshow(array, cmap=DOT_CMAP, vmin=0, vmax=1,
              interpolation='nearest', aspect='auto')
    if title:
        ax.set_title(title, fontsize=10)
    _label_axes(ax, array, s1, s2)
    return ax


def _draw_path(ax, path):
    xs = [j for i, j in path]
    ys = [i for i, j in path]
    ax.plot(xs, ys, color=ORANGE, linewidth=1.6, zorder=3)
    ax.scatter([xs[0], xs[-1]], [ys[0], ys[-1]], color=ORANGE, s=18, zorder=4)


def plot_dotplot(dotplot, s1=None, s2=None, title='Dot-plot'):
    """Dibuja el dot-plot (3.1b): s1 en el eje vertical, s2 en el horizontal."""
    fig, ax = plt.subplots(figsize=(7, 7))
    _draw_dotplot(ax, dotplot, s1, s2, title)
    fig.tight_layout()
    plt.show()


def plot_dotplot_with_path(dotplot, path, s1=None, s2=None,
                           title='Dot-plot + camino'):
    """El dot-plot con el camino superpuesto (3.1e)."""
    fig, ax = plt.subplots(figsize=(7, 7))
    _draw_dotplot(ax, dotplot, s1, s2, title)
    _draw_path(ax, path)
    fig.tight_layout()
    plt.show()


def plot_dotplot_grid(dotplots, title='Comparacion de filtros', paths=None):
    """Varios dot-plots en una grilla, para comparar filtros o umbrales (3.1c).

    dotplots: dict {etiqueta: matriz}. paths: dict opcional {etiqueta: camino}.
    """
    labels = list(dotplots)
    columns = 2 if len(labels) > 1 else 1
    rows = (len(labels) + columns - 1) // columns

    fig, axes = plt.subplots(rows, columns, figsize=(6 * columns, 5.5 * rows))
    axes = np.atleast_1d(axes).ravel()

    for ax, label in zip(axes, labels):
        _draw_dotplot(ax, dotplots[label], title=label)
        if paths and label in paths:
            _draw_path(ax, paths[label])
    for ax in axes[len(labels):]:
        ax.axis('off')

    fig.suptitle(title)
    fig.tight_layout()
    plt.show()


def plot_comparison_matrix(comparison, s1=None, s2=None,
                           title='Matriz de comparacion', path=None):
    """Mapa de calor de la matriz de puntajes de proteinas (3.2a / 3.2d).

    Escala divergente centrada en 0: azul los pares con puntaje positivo
    (sustituciones conservativas), rojo los negativos.
    """
    array = np.asarray(comparison, dtype=float)
    limit = max(abs(array.min()), abs(array.max())) or 1.0

    fig, ax = plt.subplots(figsize=(7.5, 7))
    image = ax.imshow(array, cmap='RdBu', vmin=-limit, vmax=limit,
                      interpolation='nearest', aspect='auto')
    ax.set_title(title, fontsize=10)
    _label_axes(ax, array, s1, s2)
    if path:
        _draw_path(ax, path)
    fig.colorbar(image, ax=ax, shrink=0.8, label='puntaje')
    fig.tight_layout()
    plt.show()
