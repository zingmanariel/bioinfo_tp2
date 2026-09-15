"""3.1f) Analice sobre el camino el efecto de utilizar secuencias con diferente
grado de identidad, diferente largo, y el uso de los filtros programados en el
punto anterior.

Para poder mover la identidad y el largo a voluntad, los pares de este
ejercicio se generan con sequences.mutate() y sequences.insert_indels() en vez
de leerse de data/: asi se sabe de antemano cuanta identidad tiene el par y se
puede ver que tan bien la recupera el camino.
"""

from alignment import align_from_path, alignment_stats
from dotplot import build_dotplot, density, filter_dotplot
from paths import find_path, step_counts
from plots import plot_dotplot_grid
from sequences import (generate_random_nt_sequence, insert_indels, mutate)

LENGTH = 300
IDENTITIES = [0.95, 0.85, 0.70, 0.50]
INDELS = 6
WINDOW = 10
THRESHOLD = 6


def measure(s1, s2, dotplot):
    """Corre el camino sobre el dot-plot y devuelve el resumen del alineamiento."""
    path = find_path(dotplot)
    row1, row2 = align_from_path(s1, s2, path)
    stats = alignment_stats(row1, row2)
    stats['gap_steps'] = step_counts(path)['V'] + step_counts(path)['H']
    return path, stats


def row(label, dotplot, stats):
    print(f'{label:<22} {density(dotplot):>7.2%} {stats["identity"]:>11.1%} '
          f'{stats["length"]:>8} {stats["gap_steps"]:>7}')


def header(title):
    print(f'\n--- {title} ---')
    print(f'{"":<22} {"dots":>7} {"identidad":>11} {"largo":>8} {"gaps":>7}')


def run():
    # 1) identidad: mismo largo, cada vez mas mutaciones
    header(f'efecto de la identidad (largo {LENGTH}, sin filtro)')
    dotplots, paths = {}, {}
    for identity in IDENTITIES:
        s1 = generate_random_nt_sequence(LENGTH)
        s2 = mutate(s1, identity)
        dotplot = build_dotplot(s1, s2)
        path, stats = measure(s1, s2, dotplot)
        row(f'identidad {identity:.0%}', dotplot, stats)
        label = f'identidad {identity:.0%}'
        dotplots[label], paths[label] = dotplot, path

    # 2) largo distinto: mismas mutaciones, mas indels
    header(f'efecto del largo distinto ({INDELS} indels)')
    for identity in IDENTITIES[:2]:
        s1 = generate_random_nt_sequence(LENGTH)
        s2 = insert_indels(mutate(s1, identity), INDELS)
        dotplot = build_dotplot(s1, s2)
        path, stats = measure(s1, s2, dotplot)
        row(f'identidad {identity:.0%} + indels', dotplot, stats)

    # 3) filtro: el mismo par, con y sin filtrar
    header(f'efecto del filtro (w = {WINDOW}, umbral = {THRESHOLD})')
    for identity in IDENTITIES:
        s1 = generate_random_nt_sequence(LENGTH)
        s2 = mutate(s1, identity)
        dotplot = build_dotplot(s1, s2)
        filtered = filter_dotplot(dotplot, WINDOW, THRESHOLD, diagonal=True)
        _, stats = measure(s1, s2, dotplot)
        row(f'identidad {identity:.0%} crudo', dotplot, stats)
        _, filtered_stats = measure(s1, s2, filtered)
        row(f'identidad {identity:.0%} filtrado', filtered, filtered_stats)

    print("""
A medida que baja la identidad la diagonal se despuebla, pero el camino sigue
saliendo derecho: siempre hay un camino, y con 25% de matches de fondo alcanza
para ir en diagonal aunque las secuencias no tengan nada que ver. Por eso la
identidad del alineamiento reproduce casi exactamente la que se pidio al
mutar: el camino no "descubre" homologia, la mide.

Los indels son los que sacan al camino de la diagonal: ahi aparecen los pasos
V/H, y la identidad baja mas de lo que corresponde a las mutaciones porque el
camino tarda unos pasos en reengancharse con la diagonal corrida.

El filtro no cambia el trazado mientras la identidad es alta (la diagonal
sobrevive entera y solo se apaga el fondo). Cuando la identidad es baja el
filtro tambien borra tramos de la diagonal verdadera, quedan pedazos sueltos y
el camino empieza a saltar de uno a otro: por eso la ultima fila de la tabla
tiene mas gaps filtrada que cruda.
""")

    plot_dotplot_grid(dotplots, title=f'Camino vs. identidad (largo {LENGTH})',
                      paths=paths)
