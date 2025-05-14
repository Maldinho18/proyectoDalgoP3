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
    # Construcción del grafo
    for mer in kmers:
        u = get_id(mer[:-1])
        v = get_id(mer[1:])
        adj[u].append((v, mer[-1]))
        out_deg[u] += 1
        in_deg[v] += 1

    # Orden lexicográfico de aristas
    for u in adj:
        adj[u] = deque(sorted(adj[u], key=lambda x: x[1]))

    # Encuentra nodo inicial
    start = None
    for u in range(len(rev_map)):
        if out_deg[u] - in_deg[u] == 1:
            start = u
            break
    if start is None:
        start = get_id(kmers[0][:-1])

    # Hierholzer para camino Euleriano
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


def full_superstring(kmers):
    """Reconstruye texto cubridor mínimo para todos los kmers."""
    # Conectividad simple (grafo no dirigido)
    id_map, rev_map = {}, []
    def get_id(s):
        if s not in id_map:
            id_map[s] = len(rev_map)
            rev_map.append(s)
        return id_map[s]

    undirected = defaultdict(list)
    for mer in kmers:
        u = get_id(mer[:-1])
        v = get_id(mer[1:])
        undirected[u].append(v)
        undirected[v].append(u)

    seen = [False] * len(rev_map)
    comp_strs = []
    # Recorrer componentes
    for u in range(len(rev_map)):
        if not seen[u]:
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
            # Filtrar kmers de esta componente
            sub = [mer for mer in kmers if id_map[mer[:-1]] in comp_nodes]
            comp_strs.append(eulerian_superstring(sub))

    # Concatenar directamente (sin merge costoso)
    # Componentes no solapan entre sí
    return ''.join(comp_strs)


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
