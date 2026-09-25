import pickle, math, glob
import numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
rows=[r for r in pickle.load(open('fast/rows.pkl','rb')) if r['N']>=600]
cm4={0:('#1b6ca8','o'),1:('#e0892b','s'),2:('#3a9d5d','^'),3:('#b0413e','D')}
fig,ax=plt.subplots(1,3,figsize=(17,5))
for k,(c,m) in cm4.items():
    R=[r for r in rows if r['n']%4==k]
    ax[0].errorbar([r['n'] for r in R],[r['mv'] for r in R],yerr=[1.5/math.sqrt(r['N']*20) for r in R],fmt=m,color=c,ms=6,label=f'n ≡ {k} (mod 4)',capsize=0,lw=1)
    ax[1].plot([r['n'] for r in R],[r['var'] for r in R],m,color=c,ms=6,label=f'n ≡ {k} (mod 4)')
    ax[2].plot([r['n'] for r in R],[r['tot'] for r in R],m,color=c,ms=7,label=f'n ≡ {k} (mod 4)')
nn=np.linspace(6,80,100)
ax[0].plot(nn,0.35*(nn/20)**-0.25,'k--',lw=1.2,label='∝ n^(−1/4)')
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('n'); ax[0].set_ylabel('amplitude of the 2×2-superblock pattern in ⟨v⟩\n(central disk, Fourier weight at (±π/2, ±π/2))')
ax[0].set_xticks([6,10,20,40,80]); ax[0].set_xticklabels(['6','10','20','40','80']); ax[0].set_yticks([.2,.3,.4,.5]); ax[0].set_yticklabels(['0.2','0.3','0.4','0.5']); ax[0].legend(frameon=False,fontsize=9)
ax[1].plot(nn,4.9+8/math.pi**2*np.log(nn/40),'k--',lw=1.2,label='slope 8/π² (≙ exponent 1/4)')
ax[1].set_xscale('log'); ax[1].set_xlabel('n'); ax[1].set_ylabel('Var h(centre)'); ax[1].set_xticks([6,10,20,40,80]); ax[1].set_xticklabels(['6','10','20','40','80']); ax[1].legend(frameon=False,fontsize=9)
ax[2].set_xlabel('n'); ax[2].set_ylabel('Σ v  (identical for every tiling)'); ax[2].set_yticks([-2,0,2]); ax[2].axhline(0,color='#bbb',lw=.8)
for a in ax: a.grid(alpha=.25)
plt.tight_layout(); plt.savefig('outputs/vorticity_scaling.png',dpi=190); plt.savefig('outputs/vorticity_scaling.pdf')
# ---- central maps for four residues
cmapv=LinearSegmentedColormap.from_list('ob',['#2f6db5','#e6e6e6','#e8743a']); nv=Normalize(-.5,.5)
fig,ax=plt.subplots(2,4,figsize=(19,10))
for j,n in enumerate([36,37,38,39]):
    D=pickle.load(open(f'fast/f_{n}_1500.pkl','rb')); x0,y0=D['x0'],D['y0']; cx,cy=D['center']
    V=D['vmean']; E=D['emean']; w=8
    img=np.full((2*w,2*w),np.nan)
    for i in range(-w,w):
        for k in range(-w,w):
            c=(cx+i,cy+k)
            if c in D['cells']: img[k+w,i+w]=V[c[0]-x0,c[1]-y0]
    ax[0][j].imshow(img,origin='lower',cmap=cmapv,norm=nv,extent=(cx-w,cx+w,cy-w,cy+w))
    for g in range(cx-w,cx+w+1,4): ax[0][j].axvline(g,color='k',lw=.4,alpha=.4)
    for g in range(cy-w,cy+w+1,4): ax[0][j].axhline(g,color='k',lw=.4,alpha=.4)
    ax[0][j].plot(cx,cy,'k+',ms=12,mew=2); ax[0][j].set_title(f'n = {n} (n≡{n%4} mod 4): ⟨v⟩',fontsize=12)
    img=np.full((2*w+1,2*w+1),np.nan)
    for i in range(-w,w+1):
        for k in range(-w,w+1):
            img[k+w,i+w]=E[cx+i-x0,cy+k-y0]
    ax[1][j].imshow(img,origin='lower',cmap=cmapv,norm=Normalize(-.7,.7),extent=(cx-w-.5,cx+w+.5,cy-w-.5,cy+w+.5))
    ax[1][j].plot(cx,cy,'k+',ms=12,mew=2); ax[1][j].set_title(f'n = {n}: ⟨ε⟩ (centre {E[cx-x0,cy-y0]:+.2f})',fontsize=12)
    for a in ax[:,j]: a.set_xticks([]); a.set_yticks([])
cb=fig.colorbar(plt.cm.ScalarMappable(norm=nv,cmap=cmapv),ax=ax[0].tolist(),shrink=.6,pad=.02); cb.set_label('⟨v⟩')
cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-.7,.7),cmap=cmapv),ax=ax[1].tolist(),shrink=.6,pad=.02); cb.set_label('⟨ε⟩ = ⟨(−1)^⌊h/4⌋⟩')
fig.suptitle('Centre of the Aztec diamond (+ = centre).  Top: ⟨v⟩ per cell (central 16×16). Bottom: ⟨ε⟩ per vertex. Same period-4 (2×2-superblock) wavevector for all n; n mod 4 sets its phase',fontsize=12)
plt.savefig('outputs/vorticity_centre_mod4.png',dpi=160); plt.savefig('outputs/vorticity_centre_mod4.pdf')
