import pickle, sys, math, cmath, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import hsv_to_rgb
from core import rect, pcen
from plot_defect_maps import outline, frame
AREF=0.45
def colour(z):
    ph=(cmath.phase(z)/(2*math.pi))%1.0; a=min(abs(z)/AREF,1.0)
    rgb=hsv_to_rgb((ph,0.85,0.95)); return tuple(1-(1-rgb)*a)
def panel(a,D,P,defect,title,arrows=True):
    X,Y=frame(D); Wd=10
    for c,z in P.items():
        px,py=pcen(c)
        if abs(px-X)>Wd or abs(py-Y)>Wd: continue
        a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=colour(z),ec='#888',lw=.25))
    if arrows:
        for c,z in P.items():
            if (c[0]+c[1])%2: continue
            px,py=pcen(c)
            if abs(px-X)>Wd-.6 or abs(py-Y)>Wd-.6 or abs(z)<.03: continue
            u=z/abs(z)*0.75
            a.arrow(px-u.real/2,py-u.imag/2,u.real,u.imag,head_width=.32,head_length=.3,length_includes_head=True,color='k',lw=.7)
    if defect=='O':
        for O in D['Os']:
            for c in O: a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc='#333',ec='k',lw=.8))
    if defect=='mono':
        A,B=D['AB']
        for c in (A,B): a.add_patch(Polygon(outline(list(rect(c))),closed=True,fc='#333',ec='k',lw=.8)); a.plot(*pcen(c),'*',color='w',ms=9)
    a.set_xlim(X-Wd,X+Wd); a.set_ylim(Y-Wd,Y+Wd); a.set_aspect('equal'); a.set_xticks([]); a.set_yticks([]); a.set_title(title,fontsize=11)
def wheel(ax):
    th=np.linspace(0,2*np.pi,361); r=np.linspace(0,1,40); T,Rr=np.meshgrid(th,r)
    C=np.zeros(T.shape+(3,))
    for i in range(T.shape[0]):
        for j in range(T.shape[1]): C[i,j]=colour(Rr[i,j]*AREF*cmath.exp(1j*T[i,j]))
    ax.pcolormesh(T,Rr,C[:-1,:-1],shading='flat'); ax.set_yticks([]); ax.set_xticks([k*np.pi/4 for k in range(8)])
    ax.set_xticklabels(['0\n(r=0)','π/4\n(r=1)','π/2\n(r=2)','3π/4\n(r=3)','π','5π/4','3π/2','7π/4'],fontsize=8)
    ax.set_title('arg Ψ = π h̄ / 4   (brightness: |Ψ|)',fontsize=9)
for pl in ('I','II'):
    M=pickle.load(open(f'exact/psi_{pl}.pkl','rb'))
    fig=plt.figure(figsize=(18,14.2)); gs=fig.add_gridspec(3,5,width_ratios=[1,1,1,1,.55])
    for j,n in enumerate(sorted(M)):
        D=M[n]
        panel(fig.add_subplot(gs[0,j]),D,D['none'][0],None,f'n = {n} (n ≡ {n%4}): no defect')
        panel(fig.add_subplot(gs[1,j]),D,D['O'][0],'O','O pair: two full vortices')
        panel(fig.add_subplot(gs[2,j]),D,D['mono'][0],'mono','monomer pair: two half vortices')
    wheel(fig.add_subplot(gs[0,4],projection='polar'))
    tx=fig.add_subplot(gs[1:,4]); tx.axis('off')
    tx.text(0,1,"Ψ(cell) = ½ Σ_corners e^{iπr(p)/4} ⟨ε(p)⟩\n\n⟨ε(p)⟩ = Re[e^{−iπr(p)/4} Ψ]\n⟨v(d)⟩ = s(d) Re[e^{−iπr(p_d)/4} Ψ]\n\nThe ⟨v⟩ pattern is the projection of\none complex field Ψ ≈ A e^{iπh̄/4}.\nAround an O tile Ψ turns by 2π\n(full vortex, pinwheel);\naround a monomer by π (half vortex):\nthe colour jumps to the opposite\nhue along a line = the string.\n\nArrows: direction of Ψ (every\nsecond cell). Defects 7 cells apart.",va='top',fontsize=10.5)
    fig.suptitle(f'Exact vorticity phase field Ψ, Aztec placement {pl} (Kasteleyn, no sampling)',fontsize=14)
    plt.tight_layout(rect=(0,0,1,.97))
    plt.savefig(f'outputs/psi_maps_{pl}.png',dpi=120); plt.savefig(f'outputs/psi_maps_{pl}.pdf'); print('ok',pl)
