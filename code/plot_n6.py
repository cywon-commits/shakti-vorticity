import random, math
exec(open('plot_stats.py').read().split('# ---------------- Figure 1')[0])
from matplotlib.lines import Line2D
def sample_panel(ax,n,seed):
    cells,sh=aztec(n); rng=random.Random(seed)
    d=[((a[0]+sh[0],a[1]+sh[1]),(b[0]+sh[0],b[1]+sh[1])) for a,b in shuffle(n,rng)]
    spin,free,ring,cf=analyse(d,cells)
    short=set(spin)|free
    for x in d:
        sqs=[p for c in x for p in rect(c)]
        for c in x:
            fc='#f3f3f3'
            if c in short: fc=cmap(norm(spin[c])) if c in spin else '#bdbdbd'
            ax.add_patch(Polygon(outline(list(rect(c))),closed=True,fc=fc,ec='white',lw=.5))
        ax.add_patch(Polygon(outline(sqs),closed=True,fill=False,ec='#555',lw=1))
    P={c:pcen(c) for c in short}
    for c in short:
        for dd in [(1,0),(0,1)]:
            nb=(c[0]+dd[0],c[1]+dd[1])
            if nb in short: ax.plot(*zip(P[c],P[nb]),color='k',lw=1.8,zorder=6)
    for a,b in zip(ring,ring[1:]+ring[:1]):
        ax.plot(*zip(P[a],P[b]),color='#1f77b4',ls='--',lw=1.6,zorder=5)
    for c in ring: ax.plot(*P[c],'o',ms=5,mfc='none',mec='k',mew=1.2,zorder=7)
    xs=[p[0] for c in cells for p in rect(c)]; ys=[p[1] for c in cells for p in rect(c)]
    ax.set_xlim(min(xs)-.5,max(xs)+1.5); ax.set_ylim(min(ys)-.5,max(ys)+1.5); ax.set_aspect('equal'); ax.axis('off')
    tsum=sum(spin[c]*chi(c) for c in ring)
    return ring,tsum
fig,ax=plt.subplots(2,4,figsize=(21,11))
for row,(n,mode) in enumerate([(5,'exact'),(6,'40000')]):
    D=load(n,mode); A,N=D['A'],D['N']
    ring,tsum=sample_panel(ax[row][0],n,11)
    Mt=np.mean([g['Mtau'] for g in D['glob']])
    lab='exact, all 32 768 tilings' if mode=='exact' else f'{N} uniform samples'
    ax[row][0].set_title(f'n = {n}: one random tiling\nring Σχσ = {tsum:+d}  (grey = short arm not linked to ring)',fontsize=11)
    shakti_map(ax[row][1],D,lambda c:A['sig'][c]/N,cmap,norm,ring,f'n = {n}:  ⟨σ⟩')
    shakti_map(ax[row][2],D,lambda c:A['tau'][c]/N,cmap,norm,ring,f'n = {n}:  ⟨τ⟩ = ⟨σχ⟩,   ⟨Σ_interior τ⟩ = {Mt:+.2f}')
    shakti_map(ax[row][3],D,lambda c:A['conn'][c]/N,seq,Normalize(0,1),ring,f'n = {n}:  P(linked to ring)   [{lab}]')
cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=ax[:,:3].ravel().tolist(),shrink=.45,orientation='horizontal',pad=.03); cb.set_label('purple −1  ·  yellow +1')
cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=seq),ax=ax[:,3].tolist(),shrink=.8,orientation='horizontal',pad=.03); cb.set_label('probability')
h=[Line2D([],[],color='k',lw=1.8,label='interior link'),Line2D([],[],color='#1f77b4',ls='--',lw=1.6,label='boundary-ring link'),
   Line2D([],[],marker='o',ls='',mfc='none',mec='k',label='frozen ring plaquette')]
fig.legend(handles=h,loc='upper center',ncol=3,frameon=False,fontsize=10,bbox_to_anchor=(.5,.995))
plt.savefig('outputs/ring_stats_n5_vs_n6.png',dpi=170); plt.savefig('outputs/ring_stats_n5_vs_n6.pdf')
# ring tau pattern printout
for n in (5,6):
    cells,sh=aztec(n); d=[((a[0]+sh[0],a[1]+sh[1]),(b[0]+sh[0],b[1]+sh[1])) for a,b in shuffle(n)]
    sp,fr,ring,_=analyse(d,cells); print(n,''.join('+' if sp[c]*chi(c)>0 else '-' for c in ring),''.join('H' if black(c) else 'V' for c in ring))
