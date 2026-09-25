import io, contextlib
from collections import Counter, defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    import otet
from otet import *
def stats(cells):
    deg={'D':[], 'O':[]}; isl={'D':[], 'O':[]}
    for tiles in enum_with_defects(cells):
        if not any(t[0]=='O' for t in tiles): continue
        owner=defaultdict(list)
        for k,t in enumerate(tiles):
            for e in tile_template(t): owner[e].append(k)
        nb=defaultdict(set); shared=Counter()
        for e,ks in owner.items():
            for a in ks:
                for b in ks:
                    if a!=b: nb[a].add(b)
            if len(set(ks))>1:
                for a in set(ks): shared[a]+=1
        for k,t in enumerate(tiles):
            deg[t[0]].append(len(nb[k])); isl[t[0]].append(shared[k])
    return {k:(sum(v)/len(v), sum(isl[k])/len(isl[k])) for k,v in deg.items()}
for ext in (False,True):
    otet.USE_EXT=ext; otet._cache.clear()
    print('O externals' if ext else 'O without externals', '6x6 board: (mean # coupled neighbour tiles, mean # constrained islands) =', stats({(i,j) for i in range(6) for j in range(6)}))
