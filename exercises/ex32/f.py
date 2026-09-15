"""3.2f) Utilice los modulos d) y e) para realizar nuevos alineamientos y compare
los resultados obtenidos con respecto a los anteriores (3.2c) y con los
obtenidos por BLAST (u otro similar).

Punto opcional del TP.
"""

from alignment import align_from_path, alignment_stats
from dotplot import (build_comparison_matrix, comparison_to_dotplot,
                     comparison_to_scored_dotplot, density)
from paths import find_path, find_path_scored, path_score
from sequences import (DISTANT_PROTEIN, data_path, load_fasta,
                       load_protein_pair)
from substitution import load_matrix

WINDOW = 7
THRESHOLDS = [0, 5, 8, 12, 20]   # 0 = sin filtrar, solo para tener la referencia
GAP_PENALTY = -4.0


def summarize(s1, s2, comparison, path):
    """Resumen de un camino: identidad del alineamiento y puntaje BLOSUM.

    El puntaje se mide siempre sobre la matriz de comparacion completa, no
    sobre el dot-plot filtrado: es la unica forma de comparar con la misma vara
    dos caminos que se calcularon con matrices distintas.
    """
    stats = alignment_stats(*align_from_path(s1, s2, path))
    stats['score'] = path_score(comparison, path)
    return stats


def compare(s1, s2, comparison, threshold):
    """El mismo par alineado con el dot-plot de 1s y 0s y con el de valores."""
    binary = comparison_to_dotplot(comparison, WINDOW, threshold)
    scored = comparison_to_scored_dotplot(comparison, WINDOW, threshold)
    return {
        'dots': density(binary),
        '1s y 0s': summarize(s1, s2, comparison, find_path(binary)),
        'valores': summarize(s1, s2, comparison,
                             find_path_scored(scored, gap_penalty=GAP_PENALTY)),
    }


def run():
    matrix = load_matrix()
    beta, delta = load_protein_pair()
    myoglobin = load_fasta(data_path(DISTANT_PROTEIN))

    pairs = {
        'beta vs. delta (parecidas)': (beta, delta),
        'beta vs. mioglobina (lejanas)': (beta, myoglobin),
    }

    print(f'w = {WINDOW}, gap = {GAP_PENALTY:.0f}   '
          f'(umbral 0 = matriz de comparacion sin filtrar)')

    for label, (s1, s2) in pairs.items():
        comparison = build_comparison_matrix(s1, s2, matrix)
        print(f'\n--- {label} ---')
        print(f'{"umbral":>6}{"dots":>8}   '
              f'{"1s y 0s (3.2c)":>28}   {"valores reales (3.2e)":>28}')
        print(f'{"":>6}{"":>8}   {"id":>8}{"gaps":>6}{"puntaje":>9}   '
              f'{"":>5}{"id":>8}{"gaps":>6}{"puntaje":>9}')
        for threshold in THRESHOLDS:
            result = compare(s1, s2, comparison, threshold)
            binary, scored = result['1s y 0s'], result['valores']
            print(f'{threshold:>6}{result["dots"]:>8.2%}   '
                  f'{binary["identity"]:>8.1%}{binary["gaps"]:>6}'
                  f'{binary["score"]:>9.0f}   '
                  f'{"":>5}{scored["identity"]:>8.1%}{scored["gaps"]:>6}'
                  f'{scored["score"]:>9.0f}')

    print("""
Lectura de la tabla:

1. El par parecido no distingue nada: con o sin filtro, con 1s y 0s o con
   valores, el alineamiento es el mismo (93% de identidad, sin gaps). Cuando la
   diagonal esta bien poblada, cualquier version del algoritmo la encuentra.

2. El par lejano si distingue, y en dos direcciones opuestas:
   - Sin filtrar (umbral 0), el camino sobre valores reales gana: usa los
     puntajes de BLOSUM62 para elegir entre diagonales parecidas y saca mejor
     puntaje que el binario.
   - Con umbral alto pierde, y bastante. La razon es que el filtro le saco
     justo lo que necesita: una casilla apagada vale 0, o sea "ni bueno ni
     malo", asi que el camino ya no tiene con que decidir. El binario es mas
     robusto al filtro porque solo le pide a la casilla que haya sobrevivido.

3. La conclusion practica es que el filtro y el tipo de dot-plot no son dos
   decisiones independientes: si se va a usar el dot-plot con valores reales
   conviene filtrar poco (umbral bajo), porque el valor de cada casilla es la
   informacion con la que trabaja el algoritmo.

Ojo con leer la columna de identidad como si fuera la nota del alineamiento: el
camino sobre valores reales maximiza puntaje BLOSUM, no cantidad de matches.
Puede preferir una sustitucion conservativa (I por V, K por R) antes que un
match aislado rodeado de sustituciones caras, y eso es exactamente lo que hace
un alineador de verdad.

Contra BLAST: para el par parecido los tres alineamientos (3.2c, 3.2e y BLAST)
coinciden practicamente residuo por residuo. Para el par lejano queda la
diferencia de fondo: BLAST es local y recorta el tramo que vale la pena,
nuestro camino es global y alinea todo, incluidos los extremos donde no hay
homologia; ademas BLAST usa gaps affines (abrir cuesta mas que extender), asi
que arma un gap largo donde nosotros ponemos varios cortos.
""")
