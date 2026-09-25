from pubstyle import *
import pickle, numpy as np, math, os
from matplotlib.colors import LogNorm
cols={32:'#9ecae1',64:'#4292c6',96:'#08519c',128:'#08306b'}
fig,ax=plt.subplots(1,4,figsize=(TWO,1.95),gridspec_kw=dict(wspace=0.5,left=0.065,right=0.99,top=0.9,bottom=0.22))
a=ax[0]; ins=a.inset_axes([.55,.55,.42,.4])
for L in (32,64,96):
    D=pickle.load(open(f'exact/corr_square_{L}.pkl','rb')); G=D['G']; p=(L//2,L//2)
    rs=[r for r in range(2,D['Rw']+1,2) if (p,(r,0)) in G]; g=[abs(G[(p,(r,0))]) for r in rs]
    a.loglog([r/L for r in rs],[x*L**.5 for x in g],'o-',color=cols[L],ms=2.5,lw=.8,label=f'$L={L}$')
    sl=[-math.log(g[i+1]/g[i])/math.log(rs[i+1]/rs[i]) for i in range(len(rs)-1)]
    ins.plot([(rs[i]*rs[i+1])**.5 for i in range(len(sl))],sl,'-',color=cols[L],lw=.8)
u=np.linspace(.02,.08,10); a.loglog(u,3.12*(u/.0625)**-.5,'k--',lw=.7)
import matplotlib.ticker as mt
a.set_yticks([2.5,3,4,5]); a.set_yticklabels(['2.5','3','4','5']); a.yaxis.set_minor_formatter(mt.NullFormatter())
a.set_xlabel('$r/L$'); a.set_ylabel('$L^{1/2}|G(r)|$'); a.legend(loc='lower left',handlelength=1)
ins.axhline(.5,color='k',ls='--',lw=.6); ins.set_xscale('log'); ins.set_ylim(0,.6); ins.tick_params(labelsize=5.5,pad=1,length=1.5); ins.set_title('local slope',fontsize=5.5,pad=1)
D=pickle.load(open('sq/sf_128.pkl','rb')); S=np.fft.fftshift(D['S']); q=np.fft.fftshift(D['q'])
from scipy.ndimage import maximum_filter
Sd=maximum_filter(S,size=7)
a=ax[1]; im=a.imshow(np.maximum(Sd,1e-2).T,origin='lower',extent=(q[0]/math.pi,q[-1]/math.pi,q[0]/math.pi,q[-1]/math.pi),norm=LogNorm(vmin=S.max()*1e-3,vmax=S.max()),cmap='magma',interpolation='nearest')
a.set_xlim(-1.3,1.3); a.set_ylim(-1.3,1.3); a.set_xlabel('$q_x/\\pi$'); a.set_ylabel('$q_y/\\pi$',labelpad=0); a.set_xticks([-1,0,1]); a.set_yticks([-1,0,1])
cax=a.inset_axes([0.86,0.58,0.05,0.38]); cb=fig.colorbar(im,cax=cax); cb.ax.tick_params(labelsize=5.5,colors='w',length=1.5,width=.4,pad=1); cb.outline.set_edgecolor('w'); cb.outline.set_linewidth(.4)
cb.ax.yaxis.set_ticks_position('left')
a=ax[2]; hts={}
for L in (32,64,128):
    D=pickle.load(open(f'sq/sf_{L}.pkl','rb')); S=D['S']; qq=D['q']; M=D['M']; Q=np.array(np.meshgrid(qq,qq,indexing='ij')); prof={}
    for q0 in [(.5,1),(-.5,-1),(.5,-1),(-.5,1),(1,.5),(-1,-.5),(1,-.5),(-1,.5)]:
        dq=np.hypot(((Q[0]/math.pi-q0[0]+1)%2-1)*math.pi,((Q[1]/math.pi-q0[1]+1)%2-1)*math.pi); k=np.rint(dq/(4*math.pi/M)).astype(int)
        for kk in range(0,M//8):
            sel=(k==kk)
            if sel.any(): prof.setdefault(kk,[]).append(S[sel].mean())
    ks=sorted(prof); kv=np.array(ks)*(4*math.pi/M); sv=np.array([np.mean(prof[x]) for x in ks]); hts[L]=sv[0]
    a.loglog(kv[1:]*L,sv[1:]/L**1.5,'o-',color=cols[L],ms=2,lw=.7,label=f'$L={L}$')
x=np.linspace(3.3,12,10); a.loglog(x,0.078*(x/3.3)**-1.5,'k--',lw=.8,label='slope $-3/2$')
a.set_xlabel('$|\\mathbf{q}-\\mathbf{q}_0|L$'); a.set_ylabel('$S/L^{3/2}$'); a.legend(loc='upper right',handlelength=1,fontsize=6)
a=ax[3]; Lx=sorted(hts); a.loglog(Lx,[hts[L] for L in Lx],'o',color=CAT[0],ms=4)
sl=np.polyfit(np.log(Lx),np.log([hts[L] for L in Lx]),1)[0]
a.loglog(Lx,[hts[Lx[0]]*(L/Lx[0])**1.5 for L in Lx],'k--',lw=.8,label=f'$L^{{3/2}}$ (fit {sl:.2f})')
a.set_xlabel('$L$'); a.set_ylabel('$S(\\mathbf{q}_0)$'); a.legend(loc='upper left'); a.set_xticks([32,64,128]); a.set_xticklabels(['32','64','128'])
import matplotlib.ticker as mt; a.xaxis.set_minor_formatter(mt.NullFormatter())
for a_,l in zip(ax,'abcd'): a_.set_title(f'({l})',loc='left',pad=2)
save(fig,'fig10_corr')
