import sys
from vort_pairs import *
_cache={}
def info_c(d):
    d=tuple(sorted(d))
    if d not in _cache: _cache[d]=info(d)
    return _cache[d]
def signed_graph(doms):
    I=[info_c(d) for d in doms]
    owner=defaultdict(list)
    for k,A in enumerate(I):
        for e,(b,kind) in A['tp'].items(): owner[e].append((k,b))
    adj=defaultdict(dict); bad=0
    for e,lst in owner.items():
        for x in range(len(lst)):
            for y in range(x+1,len(lst)):
                (i,bi),(j,bj)=lst[x],lst[y]
                if i==j: continue
                J=bi*bj
                if j in adj[i] and adj[i][j]!=J: bad+=1
                adj[i][j]=J; adj[j][i]=J
    return I,adj,bad,owner
def frustration(nn,adj):
    v={}; fr=0; comps=0
    for s in range(nn):
        if s in v: continue
        comps+=1; v[s]=1; Q=[s]
        for u in Q:
            for w,J in adj[u].items():
                if w not in v: v[w]=v[u]*J; Q.append(w)
                elif v[w]!=v[u]*J: fr+=1
    return v,fr//2,comps
if __name__=='__main__':
    from stats_run import enum
    for n in range(1,int(sys.argv[1])+1):
        it,cells=enum(n); C=Counter()
        for doms in it:
            I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
            C[('pair-conflict' if bad else '')+('frustrated' if fr else 'ok'),comps]+=1
        print(n,dict(C))
