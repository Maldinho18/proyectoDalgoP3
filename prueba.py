import sys
from collections import defaultdict


def problemaP3(n, coordenadas):

    puntos = [(coordenadas[i], coordenadas[i + 1]) for i in range(0, 2 * n, 2)]

    
    filas = defaultdict(list)   
    columnas = defaultdict(list)  

    for x, y in puntos:
        filas[y].append(x)
        columnas[x].append(y)


    for xs in filas.values():
        xs.sort()
    for ys in columnas.values():
        ys.sort()

    llavesCol = sorted(columnas.keys())
    minY = {x: columnas[x][0] for x in llavesCol}
    maxY = {x: columnas[x][-1] for x in llavesCol}

    
    verticales = []
    for x in llavesCol:
        verticales.append((x, minY[x], x, maxY[x]))

   
    horizontales = []
    for y in sorted(filas.keys()):
        xs = filas[y]   
        inicio = xs[0]
        actual = xs[0]

        for siguiente in xs[1:]:
            bloqueado = False
            for xc in llavesCol:
                
                
                if actual < xc < siguiente and (minY[xc] <= y <= maxY[xc]):
                    bloqueado = True
                    break

            if bloqueado:
                horizontales.append((inicio, y, actual, y))
                inicio = siguiente

            actual = siguiente

        
        horizontales.append((inicio, y, actual, y))


    partes = [str(len(horizontales))]
    for (x1, y1, x2, y2) in horizontales:
        partes.extend([str(x1), str(y1), str(x2), str(y2)])
    sys.stdout.write(" ".join(partes) + "\n")

    partes = [str(len(verticales))]
    for (x1, y1, x2, y2) in verticales:
        partes.extend([str(x1), str(y1), str(x2), str(y2)])
    sys.stdout.write(" ".join(partes) + "\n")

    sys.stdout.flush()


def main():
    
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
    
        problemaP3(n, coords)


if __name__ == "__main__":
    main()
