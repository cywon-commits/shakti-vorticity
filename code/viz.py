import random, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, FancyBboxPatch
from matplotlib.colors import to_rgb
exec(open('embed.py').read().split("print('inconsistencies'")[0])
inv0={c:k for k,c in coord.items()}
inv={(-c[0]-1,-c[1]-1):k for c,k in inv0.items()}   # 180° relabel so Shakti ≈ chessboard rotated 45°          # chess cell -> Shakti plaquette (a,b,k)

n=5
cells={(i,j) for i in range(-n,n) for j in range(-n,n) if abs(i+.5)+abs(j+.5)<=n}
assert len(cells)==2*n*(n+1)
# start: all-horizontal tiling, then random 2x2 flips
dom={}
for (i,j) in sorted(cells):
    if (i,j) in dom: continue
    dom[(i,j)]=(i+1,j); dom[(i+1,j)]=(i,j)
random.seed(7)
cl=sorted(cells)
for _ in range(20000):
    i,j=random.choice(cl)
    sq=[(i,j),(i+1,j),(i,j+1),(i+1,j+1)]
    if not all(c in cells for c in sq): continue
    if dom[(i,j)]==(i+1,j) and dom[(i,j+1)]==(i+1,j+1):
        dom[(i,j)]=(i,j+1);dom[(i,j+1)]=(i,j);dom[(i+1,j)]=(i+1,j+1);dom[(i+1,j+1)]=(i+1,j)
    elif dom[(i,j)]==(i,j+1) and dom[(i+1,j)]==(i+1,j+1):
        dom[(i,j)]=(i+1,j);dom[(i+1,j)]=(i,j);dom[(i,j+1)]=(i+1,j+1);dom[(i+1,j+1)]=(i,j+1)
doms=sorted({tuple(sorted((c,dom[c]))) for c in cells})
assert len(doms)==n*(n+1)
black=lambda c:(c[0]+c[1])%2==0
def dtype(d):
    b=d[0] if black(d[0]) else d[1]; w=d[1] if b==d[0] else d[0]
    return {(1,0):'E',(-1,0):'W',(0,1):'N',(0,-1):'S'}[(w[0]-b[0],w[1]-b[1])]
col={'E':'#e07b39','W':'#3a86c8','N':'#4caf6e','S':'#c8506a'}
def shade(c,black_cell):
    r=to_rgb(c); f=0.72 if black_cell else 1.0
    return tuple(x*f+(1-f)*(0 if black_cell else 1) for x in r) if black_cell else tuple(0.55*x+0.45 for x in r)

def outline(sqs):
    """boundary polygon of a union of unit squares (simply connected)"""
    E=set()
    for (x,y) in sqs:
        for e in [((x,y),(x+1,y)),((x+1,y),(x+1,y+1)),((x+1,y+1),(x,y+1)),((x,y+1),(x,y))]:
            if (e[1],e[0]) in E: E.remove((e[1],e[0]))
            else: E.add(e)
    nxt={a:b for a,b in E}; start=next(iter(nxt)); poly=[start]; p=nxt[start]
    while p!=start: poly.append(p); p=nxt[p]
    return poly

fig,(A,B)=plt.subplots(1,2,figsize=(15,7.6))
# ---------- chessboard / Aztec diamond ----------
for k,d in enumerate(doms,1):
    t=dtype(d)
    for c in d:
        A.add_patch(Polygon(outline([c]),closed=True,fc=shade(col[t],black(c)),ec='white',lw=0.6))
    A.add_patch(Polygon(outline(list(d)),closed=True,fill=False,ec='k',lw=2))
    cx=sum(c[0] for c in d)/2+.5; cy=sum(c[1] for c in d)/2+.5
    A.text(cx,cy,str(k),ha='center',va='center',fontsize=8,color='k')
A.set_xlim(-n-.5,n+.5);A.set_ylim(-n-.5,n+.5);A.set_aspect('equal');A.axis('off')
A.set_title(f'Aztec diamond, n = {n}: domino tiling ({len(doms)} dominoes)',fontsize=12)

# ---------- Shakti lattice / L-tetrominoes ----------
used=[inv[c] for c in cells]
for k,d in enumerate(doms,1):
    t=dtype(d); sqs=[]
    for c in d:
        r=rs[inv[c]]; sqs+=list(r)
        B.add_patch(Polygon(outline(list(r)),closed=True,fc=shade(col[t],black(c)),ec='none'))
    B.add_patch(Polygon(outline(sqs),closed=True,fill=False,ec='k',lw=2,zorder=3))
    el=next(p for p in sqs if (lambda nb: len(nb)==2 and nb[0][0]*nb[1][0]+nb[0][1]*nb[1][1]==0)([(q[0]-p[0],q[1]-p[1]) for q in sqs if abs(p[0]-q[0])+abs(p[1]-q[1])==1]))
    B.text(el[0]+.5,el[1]+.5,str(k),ha='center',va='center',fontsize=8,zorder=5)
# islands: every unit edge bounding a used plaquette, except the plaquette's internal edge
edges=set()
for key in used:
    for (x,y) in rs[key]:
        for e in [((x,y),(x+1,y)),((x,y+1),(x+1,y+1)),((x,y),(x,y+1)),((x+1,y),(x+1,y+1))]:
            edges.add(e)
for key in used:
    r=sorted(rs[key]); (x0,y0),(x1,y1)=r
    if y0==y1: edges.discard(((x1,y0),(x1,y0+1)))       # horizontal plaquette: remove middle vertical island
    else: edges.discard(((x0,y1),(x0+1,y1)))            # vertical plaquette: remove middle horizontal island
g=0.14; w=0.13
for (p,q) in edges:
    if p[1]==q[1]: B.add_patch(FancyBboxPatch((p[0]+g,p[1]-w/2),1-2*g,w,boxstyle='round,pad=0,rounding_size=0.06',fc='#555',ec='none',zorder=4))
    else: B.add_patch(FancyBboxPatch((p[0]-w/2,p[1]+g),w,1-2*g,boxstyle='round,pad=0,rounding_size=0.06',fc='#555',ec='none',zorder=4))
xs=[p[0] for k in used for p in rs[k]]; ys=[p[1] for k in used for p in rs[k]]
B.set_xlim(min(xs)-.6,max(xs)+1.6);B.set_ylim(min(ys)-.6,max(ys)+1.6);B.set_aspect('equal');B.axis('off')
B.set_title('Shakti lattice: corresponding L-tetromino tiling\n(grey bars = nanoislands; each L = two perpendicular plaquettes)',fontsize=12)

from matplotlib.patches import Patch
h=[Patch(fc=col[t],ec='k',label=f'{t}: white cell {s} black cell') for t,s in [('E','right of'),('W','left of'),('N','above'),('S','below')]]
h+= [Patch(fc='#777',label='dark shade = black cell  ↔  horizontal plaquette'),Patch(fc='#ddd',label='light shade = white cell  ↔  vertical plaquette')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=10)
plt.tight_layout(rect=(0,0.07,1,0.96))
import os; os.makedirs('outputs',exist_ok=True)
for ext in ('png','pdf'): plt.savefig(f'outputs/aztec_n5_shakti_L.{ext}',dpi=200)
# checks
allsq=[p for k in used for p in rs[k]]; assert len(allsq)==len(set(allsq))==4*len(doms)
for d in doms: assert isL(rs[inv[d[0]]],rs[inv[d[1]]])
print('black cells -> H plaquettes:',all((rs[inv[c]].__iter__().__next__()) and ((inv[c][0]+inv[c][1])%2==0)==black(c) for c in cells))
