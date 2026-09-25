import sys, random, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from energy import *
def outer_islands(doms, spin):
    # islands owned (as perimeter) by exactly one tetromino = on the outer boundary of the sample
    cnt=Counter()
    I=[info_c(d) for d in doms]
    for A in I:
        for e,(b,k) in A['tp'].items():
            if k=='per': cnt[e]+=1
    return [e for e,c in cnt.items() if c==1]
def cost(qv,zv,mode):
    if mode=='q2': return qv*qv
    return (abs(qv)-(zv%2))//2
def relaxed(doms,mode):
    spin,v,I=island_config(doms); q,z=charges(spin)
    outer=outer_islands(doms,spin)
    base=sum(cost(q[p],z[p],mode) for p in z)
    # fixed charge excluding outer islands
    qf=Counter(q)
    for (a,b) in outer:
        s=spin[(a,b)]; qf[b]-=s; qf[a]+=s
    # order the outer islands into cycles
    inc=defaultdict(list)
    for e in outer: inc[e[0]].append(e); inc[e[1]].append(e)
    assert all(len(l)==2 for l in inc.values())
    seen=set(); best_total=0
    affected=set(inc)
    rest=sum(cost(q[p],z[p],mode) for p in z if p not in affected)
    for e0 in outer:
        if e0 in seen: continue
        # walk the cycle
        cyc=[e0]; seen.add(e0); p=e0[1]
        while True:
            nxt=[e for e in inc[p] if e not in seen]
            if not nxt: break
            e=nxt[0]; cyc.append(e); seen.add(e); p=e[1] if e[0]==p else e[0]
        # vertices between consecutive edges
        m=len(cyc)
        def shared(e,f): return (set(e)&set(f)).pop()
        verts=[shared(cyc[i],cyc[(i+1)%m]) for i in range(m)]
        def contrib(e,s,p):   # charge contribution of edge e with orientation s (+1: e[0]->e[1]) at vertex p
            return s if p==e[1] else -s
        INF=10**9; best=INF
        for s0 in (1,-1):
            dp={s0:0}
            for i in range(1,m):
                nd={}
                for sp,c in dp.items():
                    for s in (1,-1):
                        p=verts[i-1]; qq=qf[p]+contrib(cyc[i-1],sp,p)+contrib(cyc[i],s,p)
                        val=c+cost(qq,z[p],mode)
                        if s not in nd or val<nd[s]: nd[s]=val
                dp=nd
            for sl,c in dp.items():
                p=verts[m-1]; qq=qf[p]+contrib(cyc[m-1],sl,p)+contrib(cyc[0],s0,p)
                best=min(best,c+cost(qq,z[p],mode))
        best_total+=best
    return base, rest+best_total
for n in (range(1,int(sys.argv[1])+1) if __name__=="__main__" else []):
    for label,shift in (('I',((n+1)%2,0)),('II',(n%2,0))):
        cells,_=aztec(n,shift); R=Counter()
        it=enum_cells(cells) if n<=4 else ([((a[0]+shift[0],a[1]+shift[1]),(b[0]+shift[0],b[1]+shift[1])) for a,b in shuffle(n,random.Random(n))] for _ in range(100))
        for doms in it:
            b1,r1=relaxed(doms,'exc'); b2,r2=relaxed(doms,'q2'); R[(b1,r1,b2,r2)]+=1
        print(f'n={n:2d} {label:2s} (excess units: template, relaxed boundary | Σq²: template, relaxed boundary) ->',dict(R))
