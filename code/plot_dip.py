import pickle, math, cmath, numpy as np, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from dip_mc import simulate
from square_mc import setup, vfield
from plot_psi import colour
from core import rect, short_cell
from plot_defect_maps import outline
L=48
snaps=simulate(L,'point',[1e9,0.8,0.6,0.5,0.45,0.4,0.3],3000,500,seed=21,tag='snap')
R={Lx:pickle.load(open(f'dip/fss_{Lx}_point.pkl','rb')) for Lx in (24,32,48,64)}
Ts=np.array([r['T'] for r in R[24]]); Tp=np.where(Ts>1e3,4.0,Ts)
fig=plt.figure(figsize=(21,12)); gs=fig.add_gridspec(2,4)
cols={24:'#c6dbef',32:'#6baed6',48:'#2171b5',64:'#08306b'}
a=fig.add_subplot(gs[0,0])
for Lx in R: a.plot(Tp,[r['C'] for r in R[Lx]],'o-',color=cols[Lx],ms=4,label=f'L = {Lx}')
a.set_xscale('log'); a.set_xlabel('T  (units D = μ₀m²/4πa³);  rightmost point: T = ∞'); a.set_ylabel('specific heat per cell'); a.legend(); a.grid(alpha=.3); a.set_title('(a) specific heat',fontsize=11)
a=fig.add_subplot(gs[0,1])
for Lx in R: a.plot(Tp,[np.mean(np.abs(r['vmap'])) for r in R[Lx]],'o-',color=cols[Lx],ms=4,label=f'L = {Lx}')
a.set_xscale('log'); a.set_xlabel('T'); a.set_ylabel('mean |⟨v⟩| over the central window'); a.grid(alpha=.3); a.legend(); a.set_title('(b) vorticity order: algebraic → long-range (|⟨v⟩| = 1)',fontsize=11)
a=fig.add_subplot(gs[0,2])
Ls=np.array(list(R)); x2=[-np.polyfit(np.log(Ls),np.log([R[Lx][i]['m2'] for Lx in R]),1)[0] for i in range(len(Ts))]
a.plot(Tp,x2,'o-',color='#08306b'); a.axhline(x2[0],color='grey',ls=':',label=f'T = ∞ estimate ({x2[0]:.2f}; exact 0.5)'); a.axhline(0,color='k',lw=.6)
a.set_xscale('log'); a.set_xlabel('T'); a.set_ylabel('2x from ⟨|Ψ̄|²⟩ ~ L^{−2x}  (L = 24–64)'); a.legend(fontsize=9); a.grid(alpha=.3); a.set_title('(c) effective exponent: unchanged at high T,\nstiffening near T_c, locked below T_c ≈ 0.45 D',fontsize=11)
# crystal chess view
a=fig.add_subplot(gs[0,3]); X,Y,sgn,Rr=setup(L,0,0); h=snaps[-1]['h']; v,Hd,Vd=vfield(h,X,Y,sgn,L)
for x in range(L):
    for y in range(L):
        if Hd[x,y]: a.add_patch(Rectangle((x+.08,y+.08),1.84,.84,fc='#e8743a' if v[x,y]>0 else '#2f6db5',ec='k',lw=.4))
        if Vd[x,y]: a.add_patch(Rectangle((x+.08,y+.08),.84,1.84,fc='#e8743a' if v[x,y]>0 else '#2f6db5',ec='k',lw=.4))
a.set_xlim(L/2-8,L/2+8); a.set_ylim(L/2-8,L/2+8); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); a.set_title('(d) low-T crystal (T = 0.3), chessboard view:\ndiagonal (π/2, π/2) order, flat height',fontsize=11)
for j,k in enumerate([0,2,4,6]):
    a=fig.add_subplot(gs[1,j]); h=snaps[k]['h']
    P0={}
    for x in range(L):
        for y in range(L): P0[(x,y)]=0.5*sum(cmath.exp(1j*math.pi*h[x+aa,y+bb]/4) for aa in (0,1) for bb in (0,1))
    for (x,y),z in P0.items():
        zz=0;ws=0
        for dx in range(-3,4):
            for dy in range(-3,4):
                q=(x+dx,y+dy)
                if q in P0: w=math.exp(-(dx*dx+dy*dy)/(2*1.3**2)); zz+=w*P0[q]; ws+=w
        a.add_patch(Polygon(outline(list(rect((x,y)))),closed=True,fc=colour(zz/ws*0.45/0.6),ec='none'))
    a.set_aspect('equal'); a.autoscale_view(); a.axis('off'); T=snaps[k]['T']; a.set_title(f'(e) phase field of one state, T = {"∞" if T>1e3 else T}',fontsize=11)
fig.suptitle('Dipolar interactions within the ground-state manifold (point dipoles, cutoff 6a): critical vorticity phase → locked crystal with long-range vorticity order',fontsize=13)
plt.tight_layout(rect=(0,0,1,.96)); plt.savefig('outputs/dipolar_ordering.png',dpi=110); plt.savefig('outputs/dipolar_ordering.pdf'); print('ok')
