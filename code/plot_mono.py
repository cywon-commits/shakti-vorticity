import pickle, random, math
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mono_tile import *
from core import rect, short_cell
OR,BL,ORl,BLl='#e8743a','#2f6db5','#f6c6a8','#aac6ea'
def tile_squares(t):
    if t[0]=='D':
        s=short_cell(t[1],t[2]); l=t[2] if s==t[1] else t[1]
        return list(rect(s)),list(rect(l))
    return list(rect(t[1]))+list(rect(t[2])),[]
def outline(ax,sqs,lw=1.3,color='k'):
    S=set(map(tuple,sqs))
    for (x,y) in S:
        for (dx,dy,seg) in [((0,-1),None,((x,y),(x+1,y))),((0,1),None,((x,y+1),(x+1,y+1))),((-1,0),None,((x,y),(x,y+1))),((1,0),None,((x+1,y),(x+1,y+1)))]:
            pass
    for (x,y) in S:
        for nb,seg in [((x,y-1),((x,y),(x+1,y))),((x,y+1),((x,y+1),(x+1,y+1))),((x-1,y),((x,y),(x,y+1))),((x+1,y),((x+1,y),(x+1,y+1)))]:
            if nb not in S: ax.plot([seg[0][0],seg[1][0]],[seg[0][1],seg[1][1]],color=color,lw=lw,solid_capstyle='round')
def cen(sqs): return (sum(p[0]+.5 for p in sqs)/len(sqs), sum(p[1]+.5 for p in sqs)/len(sqs))
def draw(ax,tiles,holes,title):
    J,conf=couplings(tiles); f,viol,u=min_frustration(len(tiles),J)
    C=[]
    for k,t in enumerate(tiles):
        s,l=tile_squares(t); v=2*u[k]-1
        if t[0]=='D':
            for p in s: ax.add_patch(Rectangle(p,1,1,fc=OR if v>0 else BL,ec='none'))
            for p in l: ax.add_patch(Rectangle(p,1,1,fc=ORl if v>0 else BLl,ec='none'))
        else:
            for p in s: ax.add_patch(Rectangle(p,1,1,fc='#6a3d9a' if v>0 else '#33a02c',ec='none'))
            c=cen(s); ax.text(c[0],c[1],'O',ha='center',va='center',color='white',fontsize=8,fontweight='bold')
        outline(ax,s+l); C.append(cen(s+l))
    for h in holes:
        sq=list(rect(h))
        for p in sq: ax.add_patch(Rectangle(p,1,1,fc='#bbbbbb',ec='none'))
        outline(ax,sq,lw=1.6); c=cen(sq); ax.plot(*c,'k*',ms=11,zorder=8)
    for (i,j) in viol:
        ax.plot([C[i][0],C[j][0]],[C[i][1],C[j][1]],color='#d00000',lw=3.5,zorder=9,solid_capstyle='round')
    xs=[p[0] for t in tiles for p in sum(tile_squares(t),[])]; ys=[p[1] for t in tiles for p in sum(tile_squares(t),[])]
    ax.set_xlim(min(xs)-.5,max(xs)+1.5); ax.set_ylim(min(ys)-.5,max(ys)+1.5); ax.set_aspect('equal'); ax.axis('off')
    ax.set_title(f'{title}\nminimum # violated couplings = {f}',fontsize=11)
L=14; board={(i,j) for i in range(L) for j in range(L)}; rng=random.Random(11)
a,b=(3,7),(10,7)
fig=plt.figure(figsize=(19,11.5)); gs=fig.add_gridspec(2,3,height_ratios=[1.35,1])
ax1=fig.add_subplot(gs[0,0]); ax2=fig.add_subplot(gs[0,1]); ax3=fig.add_subplot(gs[0,2])
tiles=random_tiling(board-{a,b},rng); draw(ax1,tiles,[a,b],'(a) two monomers (uncovered cells, ★)')
Oa=tuple(sorted((a,partner(a)))); bb=(9,8); Ob=tuple(sorted((bb,partner(bb))))
tiles=random_tiling(board-set(Oa)-set(Ob),rng)+[('O',)+Oa,('O',)+Ob]; draw(ax2,tiles,[],'(b) two O tetrominoes (charge ±2 at the centre)')
a2=(6,6); tiles=random_tiling(board-{a2,(a2[0]+1,a2[1])},rng); 
tiles=random_tiling(board-{(3,3)}-{(3,10)},rng) if False else tiles
# single monomer near the edge + partner at the edge: string ends on the boundary
tiles=random_tiling(board-{(6,6),(0,13)},rng); draw(ax3,tiles,[(6,6),(0,13)],'(c) monomer near the centre + monomer at the edge:\nthe string ends on the free boundary')
# energies
ax4=fig.add_subplot(gs[1,0]); ax5=fig.add_subplot(gs[1,1]); ax6=fig.add_subplot(gs[1,2])
R=pickle.load(open('mono/tile_scan.pkl','rb'))
import numpy as np
ds=sorted({d for d,_ in R['mono']})
for d,fs in R['mono']: ax4.plot([d]*len(fs),fs,'o',color='#d00000',alpha=.25,ms=5)
mean=[np.mean(sum([fs for d2,fs in R['mono'] if d2==d],[])) for d in ds]
ax4.plot(ds,mean,'-o',color='#d00000',label='monomer pair (mean over tilings)')
for dx,fs in R['O']: ax4.plot([dx+1]*len(fs),fs,'s',color='#6a3d9a',ms=7)
ax4.plot([],[],'s',color='#6a3d9a',label='O pair')
ax4.set_xlabel('separation (chessboard steps)'); ax4.set_ylabel('min # violated couplings'); ax4.legend(); ax4.set_title('(d) string length (22×22 board, random tilings)',fontsize=11); ax4.grid(alpha=.3)
for D,c in [(2,'#9ecae1'),(3,'#3182bd'),(5,'#08519c')]:
    try: res=pickle.load(open(f'mono/path24_D{D}.pkl','rb'))
    except Exception: continue
    keep=[];rm=0
    for r,e,x in res:
        if r>rm+1e-9 and r<13.3: keep.append((r,e)); rm=r
        else: break
    ax5.plot([k[0] for k in keep],[k[1] for k in keep],'-o',color=c,label=f'monomer pair, Δ = {D}')
for typ,c,m in [('z4q2','#6a3d9a','s'),('z2q','#33a02c','^')]:
    try: d=pickle.load(open(f'mono/pair12_{typ}.pkl','rb'))
    except Exception: continue
    lab={'z4q2':'z=4 charge ±2 pair (O centre)','z2q':'z=2 charge ±2 pair','z3q3':'z=3 charge ±3 pair'}[typ]
    ax5.plot([o[0] for o in d['out']],[o[1] for o in d['out']],m+'--',color=c,label=lab+', Δ = 3')
ax5.set_xlabel('separation r (Shakti lattice units)'); ax5.set_ylabel('excess energy (units of J_unhappy)'); ax5.legend(fontsize=8); ax5.grid(alpha=.3)
ax5.set_title('(e) Shakti vertex model: confined vs free',fontsize=11)
ax6.axis('off')
ax6.text(0,1,"Height winding around a point defect\n\n• monomer (one uncovered cell): Δh = ±4\n   ε = (−1)^⌊h/4⌋ changes sign once around it\n   → half vortex of θ = πh/4, ends a string of\n      violated couplings (energy ∝ length)\n\n• O tetromino: Δh = ±8\n   ε returns to itself → full vortex, no string\n   (finite energy, free)\n\nConsequence: the vorticity order can only be\nscreened by full vortices (O-type, charge ±2);\nhalf vortices are linearly confined in pairs or\npushed to the boundary.",va='top',fontsize=11.5,family='DejaVu Sans')
fig.text(.5,.005,'Top: orange/blue = v = +1/−1 (dark: strong arm, light: weak arm); purple/green O = O tile in its two states; red = violated couplings (string).\n(e) monomer pair: unhappy-vertex positions fixed by the tiling, spins optimised by ILP (24×24 torus). Charge pairs: fully free optimisation (12×12 torus). Δ = cost of type-II z=4 and charged vertices.',ha='center',fontsize=10.5)
plt.tight_layout(rect=(0,.04,1,1))
plt.savefig('outputs/monomer_strings.png',dpi=140); plt.savefig('outputs/monomer_strings.pdf')
