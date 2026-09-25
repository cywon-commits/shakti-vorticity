import sys
from vheight import *
n=int(sys.argv[1]); it,cells=enum(n)
tab=defaultdict(Counter)
for doms in it:
    I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
    spin,free,ring,_=analyse(doms,cells)
    sc={short_cell(*doms[k]):k for k in range(len(doms))}; ref=v[sc[ring[0]]]
    h=height(doms,cells)
    for k,(a,b) in enumerate(doms):
        d=tuple(sorted((a,b))); p=d[0]      # lower-left cell corner = vertex (x,y) of first cell
        hv=h[p]; q=(hv - (hv%4))//4
        s=v[k]*ref*(-1)**q
        typ='H' if d[1][0]!=d[0][0] else 'V'
        tab[(typ, p[0]%4, p[1]%4, hv%4)][s]+=1
bad=[k for k,c in tab.items() if len(c)>1]
print('n',n,'classes',len(tab),'ambiguous',len(bad))
for k in sorted(tab): print(k, dict(tab[k]))
