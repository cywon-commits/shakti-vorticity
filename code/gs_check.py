"""Is the set of 'tetromino-template' island configurations exactly the Shakti vertex-model ground-state manifold?
Periodic Shakti torus of P x P lattice units (P multiple of 4).  Vertex energies:
 z=4: type I (2-in-2-out, ins opposite) 0 ; anything else forbidden (large)
 z=3: |q|=1 with the collinear pair head-to-tail: 0 (type I');  |q|=1 collinear head-to-head/tail-to-tail: 1 (type II');  |q|=3 forbidden
 z=2: head-to-tail 0 ; otherwise forbidden."""
import sys, itertools
from collections import Counter, defaultdict
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from vort import is_island, key, template_for, RECTS
P=int(sys.argv[1]); FLIP=len(sys.argv)>2
def red(p): return (p[0]%P,p[1]%P)
def redE(e): return key((red(e[0]),red(e[1]))) if True else None
# islands on the torus
isl=[]
for x in range(P):
    for y in range(P):
        for e in (((x,y),(x+1,y)),((x,y),(x,y+1))):
            if is_island(e): isl.append(e)
idx={}
for k,e in enumerate(isl): idx[key((red(e[0]),red(e[1])))]=k
# careful: an island crossing the torus edge (x+1 = P) reduces to ((P-1,y),(0,y)); orientation: store as (start,end) in the plane
def island_id(e):
    return idx[key((red(e[0]),red(e[1])))]
# vertex incidences: for vertex p, list of (island id, +1 if island's canonical start->end points INTO p)
inc=defaultdict(list)
for k,(a,b) in enumerate(isl):
    # canonical orientation: a -> b (plane coords); s=+1 means arrow a->b
    inc[red(a)].append((k,'start',a,b)); inc[red(b)].append((k,'end',a,b))
def vertex_energy(p,s):
    L=inc[p]; z=len(L)
    ins=[]; dirs=[]
    for k,role,a,b in L:
        into = (s[k]==1 and role=='end') or (s[k]==-1 and role=='start')
        ins.append(into)
        d=(b[0]-a[0],b[1]-a[1]); dirs.append(d)   # axis of the island
    q=sum(1 if i else -1 for i in ins)
    if z==2:
        return 0 if q==0 else None
    if z==4:
        if q!=0: return None
        inn=[dirs[i] for i in range(4) if ins[i]]
        return 0 if inn[0]==inn[1] else None    # the two 'in' islands collinear (same axis)
    if z==3:
        if abs(q)!=1: return None
        # collinear pair = the two islands with the same axis
        ax=Counter(dirs); a=[i for i in range(3) if ax[dirs[i]]==2]
        return (1 if ins[a[0]]!=ins[a[1]] else 0) if FLIP else (0 if ins[a[0]]!=ins[a[1]] else 1)
    return 0
verts=list(inc)
# backtracking over islands with vertex checks when all incident islands assigned
order=sorted(range(len(isl)),key=lambda k:(isl[k][0][1],isl[k][0][0]))
pos={k:i for i,k in enumerate(order)}
last={p:max(pos[k] for k,_,_,_ in inc[p]) for p in verts}
check_at=defaultdict(list)
for p in verts: check_at[last[p]].append(p)
best=[10**9]; sols=[]
s=[0]*len(isl)
def rec(i,E):
    if E>best[0]: return
    if i==len(order):
        if E<best[0]: best[0]=E; sols.clear()
        sols.append(tuple(s)); return
    k=order[i]
    for val in (1,-1):
        s[k]=val; ok=True; e2=E
        for p in check_at[i]:
            ev=vertex_energy(p,s)
            if ev is None: ok=False; break
            e2+=ev
        if ok: rec(i+1,e2)
    s[k]=0
rec(0,0)
nz3=sum(1 for p in verts if len(inc[p])==3)
print(f'torus {P}x{P}: islands {len(isl)}, z3 vertices {nz3}, ground-state energy {best[0]} (= {best[0]/nz3:.2f} per z3), degeneracy {len(sols)}')
# ---- tetromino-template configurations on the same torus
# plaquettes of the torus
plaqs=[]
for a in range(P//2):
    for b in range(P//2):
        for k in (0,1): plaqs.append(frozenset(RECTS(a,b)[k]))
def redr(r): return frozenset(red(p) for p in r)
pl_red=[redr(r) for r in plaqs]
def lshape(r1,r2):
    u=r1|r2
    if len(u)!=4: return False
    xs={p[0] for p in u}; ys={p[1] for p in u}
    if sorted((len(xs),len(ys)))!=[2,3]: return False
    return sum(1 for p in u for q in u if abs(p[0]-q[0])+abs(p[1]-q[1])==1)//2==3
# candidate L's: plaquette i (plane rep) with a plane translate of plaquette j
cands=[]
for i,ri in enumerate(plaqs):
    for j,rj in enumerate(plaqs):
        if i>=j: continue
        for tx in (-P,0,P):
            for ty in (-P,0,P):
                rt=frozenset((p[0]+tx,p[1]+ty) for p in rj)
                if lshape(ri,rt): cands.append((i,j,ri,rt))
def elbow(sqs):
    S=set(sqs)
    for p in S:
        nb=[(q[0]-p[0],q[1]-p[1]) for q in S if abs(p[0]-q[0])+abs(p[1]-q[1])==1]
        if len(nb)==2 and nb[0][0]*nb[1][0]+nb[0][1]*nb[1][1]==0: return p
def tmpl(ri,rj):
    e=elbow(ri|rj); short,long_=(ri,rj) if e in ri else (rj,ri)
    tp,fr=template_for(sorted(short),sorted(long_))
    out={}
    for (a,b),(sg,kind) in tp.items():
        if kind=='ext': continue
        k=island_id((a,b)); A,B=isl[k]
        same = (red(a)==red(A) and red(b)==red(B))
        out[k]=sg if same else -sg    # orientation relative to the island's canonical a->b
    return out
# enumerate perfect matchings of plaquettes using candidate L's, then vorticities
byp=defaultdict(list)
for c in cands: byp[c[0]].append(c); byp[c[1]].append(c)
tsols=set(); ntil=0
def match(used,chosen):
    global ntil
    free=next((i for i in range(len(plaqs)) if i not in used),None)
    if free is None:
        ntil+=1
        T=[tmpl(c[2],c[3]) for c in chosen]
        for vs in itertools.product((1,-1),repeat=len(T)):
            sp={}; ok=True
            for t,v in zip(T,vs):
                for k,sg in t.items():
                    val=v*sg
                    if sp.get(k,val)!=val: ok=False; break
                    sp[k]=val
                if not ok: break
            if ok and len(sp)==len(isl): tsols.add(tuple(sp[k] for k in range(len(isl))))
        return
    for c in byp[free]:
        o=c[1] if c[0]==free else c[0]
        if o in used: continue
        match(used|{free,o},chosen+[c])
match(frozenset(),[])
S=set(sols)
print(f'tetromino tilings on torus: {ntil}; template configurations: {len(tsols)}; all are ground states: {tsols<=S}; ground states not produced by templates: {len(S-tsols)}')
