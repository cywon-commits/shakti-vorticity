import sys, random
from vort_tiling import *
def vort(doms,cells):
    I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
    spin,free,ring,_=analyse(doms,cells)
    sc={short_cell(*doms[k]):k for k in range(len(doms))}; ref=v[sc[ring[0]]]
    return {tuple(sorted(d)):v[k]*ref for k,d in enumerate(doms)}
n=int(sys.argv[1]); cells,sh=aztec(n); rng=random.Random(5)
C=Counter(); ex=None
for t in range(200):
    doms=[tuple(sorted(((a[0]+sh[0],a[1]+sh[1]),(b[0]+sh[0],b[1]+sh[1])))) for a,b in shuffle(n,rng)]
    D=set(doms)
    # find flippable blocks
    fl=[]
    for (a,b) in doms:
        if b==(a[0]+1,a[1]) and tuple(sorted(((a[0],a[1]+1),(a[0]+1,a[1]+1)))) in D: fl.append(('h',a))
        if b==(a[0],a[1]+1) and tuple(sorted(((a[0]+1,a[1]),(a[0]+1,a[1]+1)))) in D: fl.append(('v',a))
    kind,a=rng.choice(fl)
    x,y=a
    if kind=='h': old=[((x,y),(x+1,y)),((x,y+1),(x+1,y+1))]; new=[((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]
    else: old=[((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]; new=[((x,y),(x+1,y)),((x,y+1),(x+1,y+1))]
    d2=[d for d in doms if d not in old]+new
    v1=vort(doms,cells); v2=vort(d2,cells)
    changed=[d for d in v1 if d in v2 and v1[d]!=v2[d]]
    C[len(changed)]+=1
    if ex is None and len(changed)>0: ex=(doms,d2,changed,old,new,v1,v2)
    assert sum(v1.values())==sum(v2.values())
print('n',n,'distribution of #other dominoes whose v flips after one 2x2 flip:',sorted(C.items()))
import pickle; pickle.dump(ex,open('flip_example.pkl','wb'))
