import pickle, numpy as np, math, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
cols={32:'#9ecae1',64:'#4292c6',96:'#2171b5',128:'#08306b'}
fig,ax=plt.subplots(1,4,figsize=(23,5.6))
# (a) exact correlation collapse
a=ax[0]; ins=a.inset_axes([.52,.52,.45,.42])
for L in (32,64,96,128):
    f=f'exact/corr_square_{L}.pkl'
    if not os.path.exists(f): continue
    D=pickle.load(open(f,'rb')); G=D['G']; p=(L//2,L//2)
    rs=[r for r in range(2,D['Rw']+1,2) if (p,(r,0)) in G]; g=[abs(G[(p,(r,0))]) for r in rs]
    a.loglog([r/L for r in rs],[x*L**.5 for x in g],'o-',color=cols[L],ms=4,label=f'L = {L}')
    sl=[-math.log(g[i+1]/g[i])/math.log(rs[i+1]/rs[i]) for i in range(len(rs)-1)]
    ins.plot([(rs[i]*rs[i+1])**.5 for i in range(len(sl))],sl,'o-',color=cols[L],ms=3)
u=np.linspace(.02,.08,10); a.loglog(u,3.12*(u/.0625)**-.5,'k--',lw=1,label='slope −1/2')
a.set_xlabel('r / L'); a.set_ylabel('L^{1/2} |G(r)|'); a.legend(fontsize=9,loc='lower left'); a.grid(alpha=.3,which='both')
a.set_title('(a) exact phase correlation G = ⟨e^{iπ(h(r)−h(0))/4}⟩\nsquare board, collapse with exponent 1/2',fontsize=11)
ins.axhline(.5,color='k',ls='--',lw=.8); ins.set_xscale('log'); ins.set_xlabel('r',fontsize=8); ins.set_ylabel('local slope',fontsize=8); ins.tick_params(labelsize=7); ins.set_ylim(0,.6)
# (b) S(q) map
Ls=[L for L in (32,64,128) if os.path.exists(f'sq/sf_{L}.pkl')]; Lb=max(Ls)
D=pickle.load(open(f'sq/sf_{Lb}.pkl','rb')); S=np.fft.fftshift(D['S']); q=np.fft.fftshift(D['q'])
a=ax[1]; im=a.imshow(np.maximum(S,1e-2).T,origin='lower',extent=(q[0]/math.pi,q[-1]/math.pi,q[0]/math.pi,q[-1]/math.pi),norm=LogNorm(vmin=S.max()*1e-3,vmax=S.max()),cmap='magma')
for (qx,qy,c) in [(x,y,'c') for x in (-1.5,-.5,.5,1.5) for y in (-1,1)]+[(x,y,'lime') for x in (-1,1) for y in (-1.5,-.5,.5,1.5)]+[(0,y,'w') for y in (-1.5,-.5,.5,1.5)]:
    a.plot(qx,qy,'o',mfc='none',mec=c,ms=9,mew=.8)
a.set_xlim(-2,2); a.set_ylim(-2,2); a.set_xlabel('q_x / π'); a.set_ylabel('q_y / π'); fig.colorbar(im,ax=a,shrink=.8)
a.set_title(f'(b) vorticity structure factor S(q), L = {Lb} (MC)\ncircles: π(±½,±1) cyan, π(±1,±½) green, π(0,±½) white',fontsize=11)
# (c) scaled radial profile around q0 = pi(1/2, -1) family: average over the 4 strongest equivalent peaks
a=ax[2]; hts={}
for L in Ls:
    D=pickle.load(open(f'sq/sf_{L}.pkl','rb')); S=D['S']; qq=D['q']; M=D['M']
    Q=np.array(np.meshgrid(qq,qq,indexing='ij'))
    prof={}
    for q0 in [(.5,1),(-.5,-1),(.5,-1),(-.5,1),(1,.5),(-1,-.5),(1,-.5),(-1,.5)]:
        dq=np.hypot(((Q[0]/math.pi-q0[0]+1)%2-1)*math.pi,((Q[1]/math.pi-q0[1]+1)%2-1)*math.pi)
        k=np.rint(dq/(4*math.pi/M)).astype(int)
        for kk in range(0,M//8):
            sel=(k==kk)
            if sel.any(): prof.setdefault(kk,[]).append(S[sel].mean())
    ks=sorted(prof); kv=np.array(ks)*(4*math.pi/M); sv=np.array([np.mean(prof[x]) for x in ks])
    hts[L]=sv[0]
    a.loglog(kv[1:]*L,sv[1:]/L**1.5,'o-',color=cols[L],ms=3,label=f'L = {L}')
x=np.linspace(3.3,12,10); a.loglog(x,0.078*(x/3.3)**-1.5,'k--',lw=1.2,label='slope −3/2')
a.set_xlabel('|q − q₀| L'); a.set_ylabel('S / L^{3/2}'); a.legend(fontsize=9); a.grid(alpha=.3,which='both'); a.set_title('(c) peak shape: S(q₀+k) = L^{3/2} F(kL),\nF ~ (kL)^{−3/2} tail',fontsize=11)
# (d) peak height vs L
a=ax[3]; Lx=sorted(hts); a.loglog(Lx,[hts[L] for L in Lx],'o',color='#08306b',ms=8)
if len(Lx)>1:
    sl=np.polyfit(np.log(Lx),np.log([hts[L] for L in Lx]),1)[0]
    a.loglog(Lx,[hts[Lx[0]]*(L/Lx[0])**1.5 for L in Lx],'k--',label=f'L^{{3/2}}  (fit slope {sl:.2f})')
a.set_xlabel('L'); a.set_ylabel('S(q₀) (orientation-averaged)'); a.legend(); a.grid(alpha=.3,which='both'); a.set_title('(d) peak height ~ L^{2−2x} = L^{3/2}',fontsize=11)
plt.tight_layout(); plt.savefig('outputs/correlations_sfactor.png',dpi=140); plt.savefig('outputs/correlations_sfactor.pdf'); print('ok',Ls,{L:round(hts[L],1) for L in hts})
