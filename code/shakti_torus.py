"""Shakti torus geometry, vertex energies, plaquettes and L-candidates (refactored from gs_check.py)."""
import itertools, io, contextlib
from collections import Counter, defaultdict
with contextlib.redirect_stdout(io.StringIO()):
    from vort import is_island, key, template_for, RECTS
class Torus:
    def __init__(s,P):
        s.P=P
        red=lambda p:(p[0]%P,p[1]%P); s.red=red
        s.isl=[]
        for x in range(P):
            for y in range(P):
                for e in (((x,y),(x+1,y)),((x,y),(x,y+1))):
                    if is_island(e): s.isl.append(e)
        s.idx={key((red(a),red(b))):k for k,(a,b) in enumerate(s.isl)}
        s.inc=defaultdict(list)
        for k,(a,b) in enumerate(s.isl):
            s.inc[red(a)].append((k,'start',a,b)); s.inc[red(b)].append((k,'end',a,b))
        s.verts=list(s.inc)
        s.plaqs=[]
        for a in range(P//2):
            for b in range(P//2):
                for kk in (0,1): s.plaqs.append(frozenset(RECTS(a,b)[kk]))
        s.horiz=[len({p[1] for p in r})==1 for r in s.plaqs]
        s.cands=[]
        for i,ri in enumerate(s.plaqs):
            for j,rj in enumerate(s.plaqs):
                if i>=j: continue
                for tx in (-P,0,P):
                    for ty in (-P,0,P):
                        rt=frozenset((p[0]+tx,p[1]+ty) for p in rj)
                        if lshape(ri,rt): s.cands.append((i,j,ri,rt))
    def island_id(s,e): return s.idx[key((s.red(e[0]),s.red(e[1])))]
    def tmpl(s,ri,rj,ext=False):
        e=elbow(ri|rj); short,long_=(ri,rj) if e in ri else (rj,ri)
        tp,fr=template_for(sorted(short),sorted(long_))
        out={}
        for (a,b),(sg,kind) in tp.items():
            if kind=='ext' and not ext: continue
            k=s.island_id((a,b)); A,B=s.isl[k]
            same=(s.red(a)==s.red(A) and s.red(b)==s.red(B))
            out[k]=sg if same else -sg
        return out
    def vstate(s,p,spins):
        """spins: dict/list island->+-1 ; returns (z, q, ins, dirs)"""
        L=s.inc[p]; ins=[]; dirs=[]
        for k,role,a,b in L:
            ins.append((spins[k]==1 and role=='end') or (spins[k]==-1 and role=='start'))
            dirs.append((b[0]-a[0],b[1]-a[1]))
        return len(L), sum(1 if i else -1 for i in ins), ins, dirs
def classify(z,q,ins,dirs):
    """vertex type label"""
    if z==2: return 'z2ok' if q==0 else 'z2q'
    if z==4:
        if abs(q)==4: return 'z4q4'
        if abs(q)==2: return 'z4q2'
        inn=[dirs[i] for i in range(4) if ins[i]]
        return 'z4I' if inn[0]==inn[1] else 'z4II'
    if z==3:
        if abs(q)==3: return 'z3q3'
        ax=Counter(dirs); a=[i for i in range(3) if ax[dirs[i]]==2]
        return 'z3unhappy' if ins[a[0]]!=ins[a[1]] else 'z3happy'   # head-to-tail collinear = unhappy (FLIP convention)
def lshape(r1,r2):
    u=r1|r2
    if len(u)!=4: return False
    xs={p[0] for p in u}; ys={p[1] for p in u}
    if sorted((len(xs),len(ys)))!=[2,3]: return False
    return sum(1 for p in u for q in u if abs(p[0]-q[0])+abs(p[1]-q[1])==1)//2==3
def elbow(sqs):
    S=set(sqs)
    for p in S:
        nb=[(q[0]-p[0],q[1]-p[1]) for q in S if abs(p[0]-q[0])+abs(p[1]-q[1])==1]
        if len(nb)==2 and nb[0][0]*nb[1][0]+nb[0][1]*nb[1][1]==0: return p
