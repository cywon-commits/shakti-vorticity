"""Exact pair weights: Z(O pair at distance r)/Z and Z(monomer pair)/Z on a square board -> scaling dimensions."""
import math, numpy as np, pickle
from exact_defect import *
from exact_mono import MonoRegion
L=int(__import__('sys').argv[1]) if len(__import__('sys').argv)>1 else 64
cells={(i,j) for i in range(L) for j in range(L)}; R0=Region(cells); ld0=R0.ld0
c=(L//2-8,L//2); c=c if black(c) else (c[0]-1,c[1])
res={'O':[],'M':[]}
for r in range(1,20):
    c2=(c[0]+r,c[1])
    if not black(c2):
        O1,O2=O_at(c),O_at(c2)
        if not set(O1)&set(O2):
            try:
                R=Region(cells-set(O1)-set(O2)); res['O'].append((r,R.ld0-ld0))
            except RuntimeError: pass
        try:
            R=MonoRegion(cells,c,c2); res['M'].append((r,R.ld0-ld0))
        except (RuntimeError,AssertionError): pass
for k in res:
    rs=np.array([a for a,b in res[k]]); lz=np.array([b for a,b in res[k]])
    sel=rs>=3
    sl=np.polyfit(np.log(rs[sel]),lz[sel],1)[0]
    print(k,'ln Z(r)/Z:',' '.join(f'{a}:{b:.3f}' for a,b in res[k]),'  slope d lnZ/d ln r = %.3f -> 2x = %.3f'%(sl,-sl))
pickle.dump(res,open(f'exact/pair_free_{L}.pkl','wb'))
