import sys, random, io, contextlib
from vheight import *
S_TAB={}
for line in open('vformula_n5.txt'):
    if line.startswith('('):
        k,v=line.split(') {'); k=eval(k+')'); val=int(float(v.split(':')[0])); S_TAB[(k[0],k[1],k[2])]=val
def r4(p): return (2+p[0]%2-p[1]%2+2*(p[0]%2)*(p[1]%2))%4
def height_norm(doms,cells):
    h=height(doms,cells); p0=next(iter(h)); c=(r4(p0)-h[p0])%4
    return {p:x+c for p,x in h.items()}
def v_formula(d,h):
    d=tuple(sorted(d)); p=d[0]; typ='H' if d[1][0]!=d[0][0] else 'V'
    assert h[p]%4==r4(p)
    return S_TAB[(typ,p[0]%4,p[1]%4)]*(-1)**(h[p]//4)
if __name__=='__main__':
    for n in [2,3,4,6,7,8,9,10,13,16]:
        cells,sh=aztec(n); rng=random.Random(n); res=Counter()
        for t in range(40 if n>5 else 60):
            doms=[((a[0]+sh[0],a[1]+sh[1]),(b[0]+sh[0],b[1]+sh[1])) for a,b in shuffle(n,rng)]
            I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
            h=height_norm(doms,cells)
            prod={v[k]*v_formula(d,h) for k,d in enumerate(doms)}
            res['match (up to global sign)' if len(prod)==1 else 'MISMATCH']+=1
        print(n,dict(res))
