"""Snapshots: phase field e^{i pi h/4} of single tilings (Glauber-equilibrated), no defect / O pair / monomer pair."""
import random, math, cmath, pickle
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from exact_mono import *
from mono_tile import random_tiling
from plot_psi import colour
from plot_defect_maps import outline
from core import rect, pcen, short_cell
def glauber(cells,doms,sweeps,rng):
    part={}
    for a,b in doms: part[a]=b; part[b]=a
    cl=sorted(cells); N=len(cl)
    for _ in range(sweeps*N):
        x,y=cl[rng.randrange(N)]
        a,b,c,d=(x,y),(x+1,y),(x,y+1),(x+1,y+1)
        if part.get(a)==b and part.get(c)==d: part[a]=c;part[c]=a;part[b]=d;part[d]=b
        elif part.get(a)==c and part.get(b)==d: part[a]=b;part[b]=a;part[c]=d;part[d]=c
    return sorted({tuple(sorted((a,b))) for a,b in part.items()})
def snap_height(doms,R):
    h=brute_height(doms,R); return h
def cell_psi(h,cells):
    P={}
    for c in cells:
        cs=[(c[0]+a,c[1]+b) for a in (0,1) for b in (0,1)]
        if all(p in h for p in cs): P[c]=0.5*sum(cmath.exp(1j*math.pi*h[p]/4) for p in cs)/1.0
    return P
rng=random.Random(7); n=22; pl='I'
cells,sh=aztec_cells(n,pl); cx,cy=sh[0],0; A=(cx-4,cy); B=(A[0]+7,cy)
cases=[('no defect',Region(cells),[]),('O pair (full vortices)',DefRegion(cells,[O_at(A),O_at(B)]),[O_at(A),O_at(B)]),('monomer pair (half vortices)',MonoRegion(cells,A,B),[])]
fig,ax=plt.subplots(1,3,figsize=(21,7.6))
for a,(name,R,Os) in zip(ax,cases):
    R.full=cells
    doms=[(t[1],t[2]) for t in random_tiling(R.cells,rng)]
    doms=glauber(R.cells,doms,1200,rng)
    h=snap_height(doms,R); P0=cell_psi(h,R.cells); P={}
    for c in P0:
        z=0; wsum=0
        for dx in range(-3,4):
            for dy in range(-3,4):
                q=(c[0]+dx,c[1]+dy)
                if q in P0:
                    w=math.exp(-(dx*dx+dy*dy)/(2*1.3**2)); z+=w*P0[q]; wsum+=w
        P[c]=z/wsum
    for c in cells:
        if c in P: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=colour(P[c]*0.45/0.6),ec='none'))
    for d in doms:
        a.add_patch(Polygon(outline(list(rect(d[0]))+list(rect(d[1]))),closed=True,fill=False,ec='#444',lw=.35))
    for O in Os:
        a.add_patch(Polygon(outline(list(rect(O[0]))+list(rect(O[1]))),closed=True,fc='#333',ec='k',lw=1))
    if 'monomer' in name:
        for c in (A,B): a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc='#333',ec='k',lw=1)); a.plot(*pcen(c),'*',color='w',ms=10)
    Cc=[pcen(c) for c in [(cx-1,cy-1),(cx,cy-1),(cx-1,cy),(cx,cy)]]; X=sum(p[0] for p in Cc)/4; Y=sum(p[1] for p in Cc)/4; Wd=13
    a.set_xlim(X-Wd,X+Wd); a.set_ylim(Y-Wd,Y+Wd); a.set_aspect('equal'); a.axis('off'); a.set_title(f'n = {n}, one tiling: {name}',fontsize=12)
fig.text(.5,.02,'Colour = phase of e^{iπh/4}, Gaussian-averaged over σ = 1.3 cells (same colour wheel as the exact Ψ maps; pale = small amplitude); outlines = L tetrominoes.\n An O tile carries a full colour wheel, a monomer half a wheel: the opposite hues meet along a line (the string).',ha='center',fontsize=11)
plt.tight_layout(rect=(0,.06,1,.95)); plt.savefig('outputs/psi_snapshots.png',dpi=130); plt.savefig('outputs/psi_snapshots.pdf'); print('ok')
