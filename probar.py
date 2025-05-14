#!/usr/bin/env python3
# Autores: Juan Sebastián, [tu nombre]
# ProblemaP3_Optimizado.py

import sys
from collections import defaultdict, Counter, deque

def eulerian_superstring(kmers):
    """Construye superstring mínimo de una componente conexa vía De Bruijn + Eulerian Path."""
    id_map, rev_map = {}, []
    def get_id(s):
        if s not in id_map:
            id_map[s] = len(rev_map)
            rev_map.append(s)
        return id_map[s]

    adj = defaultdict(deque)
    out_deg = Counter()
    in_deg = Counter()

    # Construcción del grafo dirigido
    for mer in kmers:
        u = get_id(mer[:-1])
        v = get_id(mer[1:])
        adj[u].append((v, mer[-1]))
        out_deg[u] += 1
        in_deg[v] += 1

    # Orden lexicográfico de aristas para elegir prefijos mínimos
    for u in adj:
        adj[u] = deque(sorted(adj[u], key=lambda x: x[1]))

    # Determinar nodo inicial para camino Euleriano
    start = None
    for u in range(len(rev_map)):
        if out_deg[u] - in_deg[u] == 1:
            start = u
            break
    if start is None:
        start = get_id(kmers[0][:-1])

    # Hierholzer: construir camino Euleriano
    stack = [start]
    char_stack = []
    path = []
    while stack:
        u = stack[-1]
        if adj[u]:
            v, c = adj[u].popleft()
            stack.append(v)
            char_stack.append(c)
        else:
            stack.pop()
            if char_stack:
                path.append(char_stack.pop())

    path.reverse()
    return rev_map[start] + ''.join(path)

def overlap(a, b):
    """Máximo solapamiento suffix(a) == prefix(b)."""
    mx = min(len(a), len(b))
    for l in range(mx, 0, -1):
        if a.endswith(b[:l]):
            return l
    return 0

def greedy_merge(strings):
    """Fusiona iterativamente el par con mayor overlap."""
    while len(strings) > 1:
        best, bi, bj = -1, 0, 1
        for i in range(len(strings)):
            for j in range(len(strings)):
                if i != j:
                    ov = overlap(strings[i], strings[j])
                    if ov > best or (ov == best and (strings[i] + strings[j][ov:]) < (strings[bi] + strings[bj][best:])):
                        best, bi, bj = ov, i, j
        merged = strings[bi] + strings[bj][best:]
        # reemplazar los dos por la cadena fusionada
        if bi < bj:
            strings.pop(bj)
            strings[bi] = merged
        else:
            strings.pop(bi)
            strings[bj] = merged
    return strings[0]

def full_superstring(kmers):
    """Reconstruye superstring mínimo para todos los kmers, manejando múltiples componentes."""
    # Mapa de ids para nodos (subcadenas de largo k-1)
    id_map, rev_map = {}, []
    def get_id(s):
        if s not in id_map:
            id_map[s] = len(rev_map)
            rev_map.append(s)
        return id_map[s]

    # Grafo no dirigido para detectar componentes
    undirected = defaultdict(list)
    for mer in kmers:
        u = get_id(mer[:-1])
        v = get_id(mer[1:])
        undirected[u].append(v)
        undirected[v].append(u)

    seen = [False] * len(rev_map)
    comp_strs = []

    # Para cada componente, decidir Euleriano o greedy según grados
    for u in range(len(rev_map)):
        if not seen[u]:
            # recolectar nodos de la componente
            stack = [u]
            comp_nodes = set()
            seen[u] = True
            while stack:
                x = stack.pop()
                comp_nodes.add(x)
                for w in undirected[x]:
                    if not seen[w]:
                        seen[w] = True
                        stack.append(w)
            # filtrar kmers de esta componente
            sub = [mer for mer in kmers if get_id(mer[:-1]) in comp_nodes]
            # calcular grados para saber si hay camino Euleriano
            out_deg = Counter()
            in_deg = Counter()
            for mer in sub:
                out_deg[mer[:-1]] += 1
                in_deg[mer[1:]] += 1
            diffs = [out_deg[p] - in_deg[p] for p in set(out_deg) | set(in_deg)]
            # Euleriano si |diff|<=1 y a lo sumo un +1 y un -1
            if diffs.count(1) <= 1 and diffs.count(-1) <= 1 and all(abs(d) <= 1 for d in diffs):
                comp_strs.append(eulerian_superstring(sub))
            else:
                comp_strs.append(greedy_merge(sub))

    # finalmente, merge de componentes de menor a mayor longitud
    comp_strs.sort(key=len)
    return greedy_merge(comp_strs)

def main():
    data = sys.stdin.buffer.read().split()
    t = int(data[0]); idx = 1
    out = []
    for _ in range(t):
        n, k = map(int, data[idx:idx+2]); idx += 2
        kmers = [b.decode() for b in data[idx:idx+n]]; idx += n
        out.append(full_superstring(kmers))
    sys.stdout.write("\n".join(out))

if __name__ == "__main__":
    main()
