from pubstyle import *
import pickle, math, numpy as np, matplotlib.ticker as mt
out=pickle.load(open('exact/scan_8_34.pkl','rb'))+pickle.load(open('exact/scan_35_58.pkl','rb'))
out=[r for r in out if r['lnerr']<1e-5]
def g(m):
    s=0.0
    for j in range(1,4*m,2): s+=min(4*m-j,3*j-4*m)*math.log(j)
    return math.exp(s-2*m*(2*m-1)*math.log(2))
fig,ax=plt.subplots(1,3,figsize=(TWO,2.15),gridspec_kw=dict(wspace=0.42,left=0.07,right=0.99,top=0.9,bottom=0.2))
for k in range(4):
    R=[r for r in out if r['n']%4==k]
    ax[0].plot([r['n'] for r in R],[r['A'] for r in R],MK[k],color=CAT[k],ms=3,label=f'$n\\equiv{k}$')
nn=np.linspace(8,60,50); ax[0].plot(nn,0.35*(nn/20)**-0.25,'k--',lw=.8,label='$\\propto n^{-1/4}$')
R=[r for r in out if r['n']%4==3]; ax[1].plot([(r['n']+1)//4 for r in R],[abs(r['eps_c']) for r in R],'o',color=CAT[0],ms=3.5,label='Kasteleyn (float)')
from eps_formula import vals
ax[1].plot(list(vals),[float(abs(v)) for v in vals.values()],'s',mfc='none',mec=CAT[1],ms=5,mew=.9,label='exact rationals')
mm=np.arange(1,16); ax[1].plot(mm,[g(m) for m in mm],'-',color='k',lw=.8,label='product formula')
mm2=np.linspace(1,15,50); C=math.exp(0.25)*2**(1/3)/1.2824271291006226**3; ax[1].plot(mm2,C*mm2**-0.25,':',color='#888',lw=1,label='$C\\,m^{-1/4}$')
E={(r['n']+1)//4:abs(r['eps_c']) for r in out if r['n']%4==3}; ms=sorted(E)
xs=[(m*(m-1))**.5 for m in ms[1:] if m-1 in E]; ys=[math.log(E[m-1]/E[m])/math.log(m/(m-1)) for m in ms[1:] if m-1 in E]
ax[2].plot(xs,ys,'o',color=CAT[0],ms=3.5,label='Kasteleyn'); mf=np.arange(2,60); ax[2].plot([(m*(m-1))**.5 for m in mf],[math.log(g(m-1)/g(m))/math.log(m/(m-1)) for m in mf],'k-',lw=.8,label='product formula')
ax[2].axhline(.25,color='#888',ls=':',lw=1)
for a in ax[:2]:
    a.set_xscale('log'); a.set_yscale('log'); a.xaxis.set_minor_formatter(mt.NullFormatter()); a.yaxis.set_minor_formatter(mt.NullFormatter())
ax[0].set_xticks([8,16,32,58]); ax[0].set_xticklabels(['8','16','32','58']); ax[0].set_yticks([.28,.32,.36,.4]); ax[0].set_yticklabels(['0.28','0.32','0.36','0.40'])
ax[0].set_xlabel('$n$'); ax[0].set_ylabel('local amplitude ($r\\leq4$)'); ax[0].legend(ncol=1,handletextpad=.3,loc='lower left')
ax[1].set_xticks([1,2,4,8,15]); ax[1].set_xticklabels(['1','2','4','8','15']); ax[1].set_yticks([.4,.5,.6,.7]); ax[1].set_yticklabels(['0.4','0.5','0.6','0.7'])
ax[1].set_xlabel('$m$'); ax[1].set_ylabel('$|\\langle\\varepsilon(\\mathrm{center})\\rangle|$'); ax[1].legend(loc='upper right',fontsize=6.5)
ax[2].set_xscale('log'); ax[2].set_xlabel('$m$'); ax[2].set_ylabel('local exponent'); ax[2].set_ylim(0.2,0.3); ax[2].legend(loc='upper right')
for a,l in zip(ax,'abc'): a.set_title(f'({l})',loc='left',pad=2)
save(fig,'fig5_exactscaling')
