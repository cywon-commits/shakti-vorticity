import numpy as np
from collections import Counter
from aztec_nb import shuffle, height, SH, SV, R4
from vform2 import height_norm, v_formula
from core import aztec
# uniformity n=3
C=Counter()
for s in range(32000):
    T=shuffle(3,s); C[T.tobytes()]+=1
v=list(C.values()); print('n=3 distinct',len(C),'chi2 %.1f (dof 63)'%sum((x-500)**2/500 for x in v))
# height/vorticity agreement with the python implementation
n=12; shx,shy=(n+1)%2,0; O=n+1; bad=0
for s in range(20):
    T=shuffle(n,s); h=height(T,n,shx,shy,R4)
    doms=[]
    for i in range(-n,n):
        for j in range(-n,n):
            t=T[i+O,j+O]
            if t in (1,2): doms.append(((i+shx,j+shy),(i+1+shx,j+shy)))
            elif t in (3,4): doms.append(((i+shx,j+shy),(i+shx,j+1+shy)))
    cells,_=aztec(n); assert {c for d in doms for c in d}==cells
    hp=height_norm(doms,cells)
    for (a,b) in doms:
        x,y=a; i,j=x-shx,y-shy
        e=1 if ((h[i+O,j+O]-(h[i+O,j+O]%4))//4)%2==0 else -1
        vn=(SH if b[1]==a[1] else SV)[x%4,y%4]*e
        if vn!=v_formula((a,b),hp): bad+=1
print('vorticity mismatches numba vs python:',bad)
# diagnose
from collections import defaultdict
for s in range(3):
    T=shuffle(n,s); h=height(T,n,shx,shy,R4); doms=[]
    for i in range(-n,n):
        for j in range(-n,n):
            t=T[i+O,j+O]
            if t in (1,2): doms.append(((i+shx,j+shy),(i+1+shx,j+shy)))
            elif t in (3,4): doms.append(((i+shx,j+shy),(i+shx,j+1+shy)))
    cells,_=aztec(n); hp=height_norm(doms,cells)
    diff=Counter()
    for (a,b) in doms:
        x,y=a; i,j=x-shx,y-shy
        diff[(h[i+O,j+O]-hp[a])]+=1
    print('sample',s,'h_numba - h_python at domino corners:',dict(diff))
