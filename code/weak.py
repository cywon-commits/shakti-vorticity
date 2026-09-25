"""Signed constraint graph on dominoes.
 node      : domino, colour sigma = colour of its short (strong) arm
 anti edge : short arms in perpendicular contact (orthogonal chess neighbours)  -> opposite colour
 anti edge : consecutive frozen-ring short arms                                   -> opposite colour
 ferro edge: long (weak) arms in perpendicular contact AND short arms not linked  -> same colour"""
import sys, math, random
from collections import Counter
from core import *
N4=[(1,0),(-1,0),(0,1),(0,-1)]
def build(doms, cells, use_ferro=True):
    spin0,free0,ring,_=analyse(doms,cells)            # ring from previous definition
    S={}; L={}
    for d in doms:
        s=short_cell(*d); l=d[0] if s==d[1] else d[1]; S[s]=s; L[l]=s   # key every domino by its short cell
    adj={s:[] for s in S}
    anti=set()
    for s in S:
        for dx,dy in N4:
            nb=(s[0]+dx,s[1]+dy)
            if nb in S: anti.add(frozenset((s,nb)))
    for a,b in zip(ring,ring[1:]+ring[:1]): anti.add(frozenset((a,b)))
    ferro=set()
    if use_ferro:
        for l,s in L.items():
            for dx,dy in N4:
                nb=(l[0]+dx,l[1]+dy)
                if nb in L and L[nb]!=s:
                    e=frozenset((s,L[nb]))
                    if e not in anti: ferro.add(e)
    for e in anti: a,b=tuple(e); adj[a].append((b,-1)); adj[b].append((a,-1))
    for e in ferro: a,b=tuple(e); adj[a].append((b,1)); adj[b].append((a,1))
    return S,adj,ring,anti,ferro
def solve(S,adj,ring):
    spin={}; frus=0
    for i,c in enumerate(ring): spin[c]=1 if i%2==0 else -1
    Q=list(ring)
    for u in Q:
        for v,sg in adj[u]:
            want=spin[u]*sg
            if v not in spin: spin[v]=want; Q.append(v)
            elif spin[v]!=want: frus+=1
    # components not reached from the ring
    comps=[]; seen=set(spin)
    for s in S:
        if s in seen: continue
        comp={s:1}; Q=[s]; seen.add(s)
        for u in Q:
            for v,sg in adj[u]:
                if v not in comp: comp[v]=comp[u]*sg; Q.append(v); seen.add(v)
                elif comp[v]!=comp[u]*sg: frus+=1
        comps.append(comp)
    return spin, comps, frus//2
if __name__=='__main__':
    n=int(sys.argv[1]); cells,sh=aztec(n); cl=sorted(cells,key=lambda c:(c[1],c[0]))
    C=Counter(); fr=[]; nf=[]
    def rec(used,ds):
        free=next((c for c in cl if c not in used),None)
        if free is None:
            S,adj,ring,anti,ferro=build(ds,cells); spin,comps,f=solve(S,adj,ring)
            C['frustrated' if f else 'ok']+=1; nf.append(len(ferro)); fr.append(len(spin)/len(S)); return
        for dd in [(1,0),(0,1)]:
            o=(free[0]+dd[0],free[1]+dd[1])
            if o in cells and o not in used: rec(used|{free,o},ds+[(free,o)])
    rec(frozenset(),[])
    print(n,dict(C),'mean #ferro %.2f'%(sum(nf)/len(nf)),'mean frac fixed by ring %.3f'%(sum(fr)/len(fr)))
