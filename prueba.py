#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from collections import defaultdict


def procesar_caso(n, coordenadas):
    """
    n      : número de focos
    coords : lista [x1, y1, x2, y2, ..., xn, yn]
    """

    puntos = [(coordenadas[i], coordenadas[i + 1]) for i in range(0, 2 * n, 2)]

    # agrupamos por filas y columnas
    filas = defaultdict(list)   # y -> lista de xs
    columnas = defaultdict(list)  # x -> lista de ys

    for x, y in puntos:
        filas[y].append(x)
        columnas[x].append(y)

    # ordenamos
    for xs in filas.values():
        xs.sort()
    for ys in columnas.values():
        ys.sort()

    llavesCol = sorted(columnas.keys())
    minY = {x: columnas[x][0] for x in llavesCol}
    maxY = {x: columnas[x][-1] for x in llavesCol}

    # 1) trayectorias verticales: una por cada x
    verticales = []
    for x in llavesCol:
        verticales.append((x, minY[x], x, maxY[x]))

    # 2) trayectorias horizontales: separamos donde una vertical cruzaría sin foco
    horizontales = []
    for y in sorted(filas.keys()):
        xs = filas[y]   # ya está ordenada
        inicio = xs[0]
        actual = xs[0]

        for siguiente in xs[1:]:
            bloqueado = False
            for xc in llavesCol:
                # ¿hay una vertical entre actual y siguiente que pase por y
                #  pero sin foco en (xc, y)? (si hubiese foco, xc estaría en xs)
                if actual < xc < siguiente and (minY[xc] <= y <= maxY[xc]):
                    bloqueado = True
                    break

            if bloqueado:
                # cerramos segmento hasta "actual"
                horizontales.append((inicio, y, actual, y))
                inicio = siguiente

            actual = siguiente

        # último tramo de la fila
        horizontales.append((inicio, y, actual, y))

    # --- imprimir este caso ---

    # línea horizontales
    partes = [str(len(horizontales))]
    for (x1, y1, x2, y2) in horizontales:
        partes.extend([str(x1), str(y1), str(x2), str(y2)])
    sys.stdout.write(" ".join(partes) + "\n")

    # línea verticales
    partes = [str(len(verticales))]
    for (x1, y1, x2, y2) in verticales:
        partes.extend([str(x1), str(y1), str(x2), str(y2)])
    sys.stdout.write(" ".join(partes) + "\n")

    # Para ver salida en consola Windows inmediatamente:
    sys.stdout.flush()


def main():
    # leemos el número de casos
    primera = sys.stdin.readline()
    if not primera:
        return
    t = int(primera.strip())

    for _ in range(t):
        linea = sys.stdin.readline()
        if not linea:
            break
        vals = list(map(int, linea.strip().split()))
        n = vals[0]
        coords = vals[1:]
        # según el enunciado, cada caso viene en UNA línea: len(coords) == 2*n
        procesar_caso(n, coords)


if __name__ == "__main__":
    main()
