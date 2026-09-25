import sys
from weak import *
def build_v(doms,cells,variant,ring_on=True):
    spin0,free0,ring,_=analyse(doms,cells)
    if not ring_on: ring=[]
    S={};L={}
    for d in doms:
        s=short_cell(*d); l=d[0] if s==d[1] else d[1]; S[s]=s; L[l]=s
    anti=set()
    for s in S:
        for dx,dy in N4:
            nb=(s[0]+dx,s[1]+dy)
            if nb in S: anti.add(frozenset((s,nb)))
    for a,b in zip(ring,ring[1:]+ring[:1]):
        if a!=b: anti.add(frozenset((a,b)))
    deg=Counter(x for e in anti for x in e)
    ferro=set()
    for l,s in L.items():
        for dx,dy in N4:
            nb=(l[0]+dx,l[1]+dy)
            if nb in L and L[nb]!=s:
                t=L[nb]; e=frozenset((s,t))
                if e in anti: continue
                if variant=='pair' or (variant=='both_isolated' and deg[s]==0 and deg[t]==0) or (variant=='one_isolated' and (deg[s]==0 or deg[t]==0)):
                    ferro.add(e)
    adj={s:[] for s in S}
    for e in anti: a,b=tuple(e); adj[a].append((b,-1)); adj[b].append((a,-1))
    for e in ferro: a,b=tuple(e); adj[a].append((b,1)); adj[b].append((a,1))
    return S,adj,ring,ferro
n=int(sys.argv[1]); cells,sh=aztec(n); cl=sorted(cells,key=lambda c:(c[1],c[0]))
R={}
def rec(used,ds):
    free=next((c for c in cl if c not in used),None)
    if free is None:
        for v in ['pair','one_isolated','both_isolated']:
            for ring_on in (True,False):
                S,adj,ring,ferro=build_v(ds,cells,v,ring_on)
                if not ring:   # need a seed for solve(): use empty ring
                    spin,comps,f=solve(S,adj,[])
                else: spin,comps,f=solve(S,adj,ring)
                k=(v,'ring' if ring_on else 'no ring'); R.setdefault(k,Counter())['frustrated' if f else 'ok']+=1
        return
    for dd in [(1,0),(0,1)]:
        o=(free[0]+dd[0],free[1]+dd[1])
        if o in cells and o not in used: rec(used|{free,o},ds+[(free,o)])
rec(frozenset(),[])
for k,v in sorted(R.items()): print(n,k,dict(v))
