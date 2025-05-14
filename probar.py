#!/usr/bin/env python3
# Autores: Juan Sebastián, [tu nombre]
# ProblemaP3.py

import sys
from collections import defaultdict, Counter, deque

def eulerian_superstring(kmers):
    """Superstring mínimo de una componente conexa vía De Bruijn + Eulerian Path."""
    id_map, rev_map = {}, []
    def get_id(s):
        if s not in id_map:
            id_map[s] = len(rev_map)
            rev_map.append(s)
        return id_map[s]

    # 1) Construcción del grafo y grados
    adj = defaultdict(list)
    out_deg = Counter()
    in_deg  = Counter()
    for mer in kmers:
        u, v = get_id(mer[:-1]), get_id(mer[1:])
        adj[u].append((v, mer[-1]))
        out_deg[u] += 1
        in_deg[v]  += 1

    # 2) Orden lexicográfico ASC y deque
    for u in adj:
        adj[u].sort(key=lambda x: x[1])
        adj[u] = deque(adj[u])

    # 3) Nodo inicial (out-in == 1) o prefijo del primer k-mer
    start = None
    for u in range(len(rev_map)):
        if out_deg[u] - in_deg[u] == 1:
            start = u
            break
    if start is None:
        start = get_id(kmers[0][:-1])

    # 4) Hierholzer
    node_stack = [start]
    char_stack = []
    path = []
    while node_stack:
        u = node_stack[-1]
        if adj[u]:
            v, c = adj[u].popleft()
            node_stack.append(v)
            char_stack.append(c)
        else:
            node_stack.pop()
            if char_stack:
                path.append(char_stack.pop())

    # 5) Reconstrucción
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
                    if ov > best:
                        best, bi, bj = ov, i, j
        merged = strings[bi] + strings[bj][best:]
        # reemplazar el par por la fusión
        if bi < bj:
            strings.pop(bj)
            strings[bi] = merged
        else:
            strings.pop(bi)
            strings[bj] = merged
    return strings[0]

def full_superstring(kmers):
    """Detecta componentes, obtiene su superstring y luego las fusiona."""
    # 1) Conectividad simple para componentes (grafo no dirigido)
    id_map, rev_map = {}, []
    def get_id_simple(s):
        if s not in id_map:
            id_map[s] = len(rev_map)
            rev_map.append(s)
        return id_map[s]

    adj = defaultdict(list)
    for mer in kmers:
        u, v = get_id_simple(mer[:-1]), get_id_simple(mer[1:])
        adj[u].append(v)
        adj[v].append(u)

    seen = [False]*len(rev_map)
    comps = []
    for u in range(len(rev_map)):
        if not seen[u]:
            stack, comp = [u], []
            while stack:
                x = stack.pop()
                if seen[x]: continue
                seen[x] = True
                comp.append(x)
                for w in adj[x]:
                    if not seen[w]:
                        stack.append(w)
            comps.append(comp)

    # 2) Superstring por componente
    comp_strs = []
    for comp in comps:
        comp_set = set(comp)
        sub = [mer for mer in kmers if get_id_simple(mer[:-1]) in comp_set]
        comp_strs.append(eulerian_superstring(sub))

    # 3) Fusionar, empezando por los más largos
    comp_strs.sort(key=len)
    return greedy_merge(comp_strs)

def main():
    data = sys.stdin.read().strip().split()
    t = int(data[0]); idx = 1
    out = []
    for _ in range(t):
        n, k = map(int, data[idx:idx+2]); idx += 2
        kmers = data[idx:idx+n]; idx += n
        out.append(full_superstring(kmers))
    sys.stdout.write("\n".join(out))

if __name__ == "__main__":
    main()
