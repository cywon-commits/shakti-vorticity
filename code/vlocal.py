import sys
from vort_tiling import *
from stats_run import enum
n=int(sys.argv[1]); it,cells=enum(n)
vals=defaultdict(Counter)
for doms in it:
    I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
    spin,free,ring,_=analyse(doms,cells)
    sc={short_cell(*doms[k]):k for k in range(len(doms))}; ref=v[sc[ring[0]]]
    for k,d in enumerate(doms): vals[tuple(sorted(d))][v[k]*ref]+=1
multi=[d for d,c in vals.items() if len(c)>1]
print('n',n,'placements used',len(vals),'placements with both signs',len(multi))
import pickle; pickle.dump({d:dict(c) for d,c in vals.items()},open(f'vloc_{n}.pkl','wb'))
