"""Finite-temperature boundary selection in block-aligned straight-edged samples (data: edge_exact.py)."""
from pubstyle import *
import pickle, numpy as np, math
NS=[11,17,23,35,47]; D={n:pickle.load(open(f'exact/edge_exact_n{n}.pkl','rb')) for n in NS}
cols=dict(zip(NS,['#c6dbef','#9ecae1','#6baed6','#2171b5','#08306b']))
fig,ax=plt.subplots(1,3,figsize=(TWO,2.15),gridspec_kw=dict(wspace=0.42,left=0.07,right=0.99,bottom=0.2,top=0.9))
# (a) low-T weight of paired edge plaquettes: <N_p> e^{J/T} -> C 2^n
import subprocess
g=[]
for n in NS:
    if n>35: g.append(np.nan); continue
    x=(n+2)*math.log(2)+8
    out=subprocess.run(['python3','edge_np.py',str(n),str(x)],capture_output=True,text=True).stdout.split()
    g.append(float(out[2])*math.exp(x))
g=np.array(g); ok=~np.isnan(g)
c=np.polyfit(np.array(NS)[ok],np.log(g[ok]),1)
a=ax[0]; a.semilogy(np.array(NS)[ok],g[ok],'o',color=CAT[0]); nn=np.linspace(8,50,10); a.semilogy(nn,2.0**(nn+2),'k--',lw=.8,label=r'$2^{\,n+2}$')
a.set_xlabel('order $n$'); a.set_ylabel(r'$\langle N_p\rangle\,e^{J/T}$ ($T\to0$)'); a.legend(loc='upper left'); a.set_title('(a) weight of one paired pair',loc='left')
print('fit slope',c[0],'ln2',math.log(2),'intercept',c[1])
for n in NS:
    o=D[n]['out']; xs=np.array(sorted(o))[::-1]; xs=xs[xs>0]
    t=n/xs
    ax[1].semilogx(t,[o[x]['paired'] for x in xs],'-o',ms=2,color=cols[n],label=f'$n={n}$')
    ax[2].semilogx(t,[o[x]['amp'] for x in xs],'-o',ms=2,color=cols[n])
ax[1].set_xlabel(r'$nT/J$'); ax[1].set_ylabel('fraction of paired edge plaquettes'); ax[1].legend(loc='upper left',fontsize=6); ax[1].set_title('(b) boundary disorder',loc='left'); ax[1].set_ylim(0,1)
ax[2].set_xlabel(r'$nT/J$'); ax[2].set_ylabel(r'central $|\Psi|/\Psi_{\rm sat}$'); ax[2].set_title('(c) imprinted order',loc='left')
for a_ in ax[1:]: a_.axvline(1/math.log(2),color='#888',ls=':',lw=.8); a_.set_xlim(0.3,60)
ax[2].text(1/math.log(2)*1.08,0.42,r'$T_{\rm on}$',fontsize=7,color='#555')
save(fig,'fig_edgeT')
