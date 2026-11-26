"""Fig. 14: finite-temperature boundary of straight-edged samples, block-aligned and through-block (data: edge_exact2.py)."""
from pubstyle import *
import pickle, numpy as np, math
from edge_half import square, EdgeModel
NS=[11,17,23,35,47]; KD=dict(aligned=4,through=2)
D={(k,n):pickle.load(open(f'exact/edge_half_{k}_n{n}.pkl','rb'))['out'] for k in KD for n in NS}
cols=dict(zip(NS,['#c6dbef','#9ecae1','#6baed6','#2171b5','#08306b']))
fig,ax=plt.subplots(1,3,figsize=(TWO,2.2),gridspec_kw=dict(wspace=0.42,left=0.07,right=0.99,bottom=0.2,top=0.9))
# (a) low-temperature weight of one excitation: <N_d> e^{J/T} -> 2^{n+c}
a=ax[0]; nn=np.array([5,7,9,11,15,19,23]); mk=dict(aligned='o',through='s')
for k,c in KD.items():
    g=[]
    for n in nn:
        cells,free,nh=square(2*n+2,k); E=EdgeModel(cells,free,nh); x=(n+c)*math.log(2)+10
        M,Mi,pf,ph=E.stats(x); g.append((len(free)-pf.sum()+ph.sum())*math.exp(x))
    a.semilogy(nn,g,mk[k],ms=3.5,mfc='white' if k=='through' else CAT[0],color=CAT[0],label=f'{"block-aligned" if k=="aligned" else "through-block"}')
    xx=np.linspace(4,24,10); a.semilogy(xx,2.0**(xx+c),'k--' if k=='aligned' else 'k:',lw=.8)
    print(k,[round(v/2**(n+c),4) for v,n in zip(g,nn)])
a.text(8.2,2.0**(8.2+4)*5,r'$2^{\,n+4}$',fontsize=7); a.text(17.5,2.0**(17.5+2)/14,r'$2^{\,n+2}$',fontsize=7)
a.set_xlabel('order $n$'); a.set_ylabel(r'$\langle N_d\rangle\,e^{J/T}$ ($T\to0$)'); a.legend(loc='upper left',fontsize=6); a.set_title('(a) boundary defects, $T\\to0$',loc='left')
for (k,n),o in D.items():
    ts=np.array(sorted(t for t in o if t<1e5)); ls='-' if k=='aligned' else '--'
    ax[1].semilogx(ts,[o[t]['Nd']/(4*(n+1)) for t in ts],ls,color=cols[n],lw=1,label=f'$n={n}$' if k=='aligned' else None)
    ax[2].semilogx(ts,[o[t]['amp'] for t in ts],ls,color=cols[n],lw=1)
ax[1].set_xlabel(r'$T/T_{\rm on}$'); ax[1].set_ylabel(r'$\langle N_d\rangle/4(n+1)$'); ax[1].legend(loc='upper left',fontsize=6); ax[1].set_title('(b) boundary disorder',loc='left'); ax[1].set_ylim(0,0.55)
ax[2].set_xlabel(r'$T/T_{\rm on}$'); ax[2].set_ylabel(r'central $|\Psi|/\Psi_{\rm sat}$'); ax[2].set_title('(c) imprinted order',loc='left'); ax[2].set_ylim(0,0.5)
for a_ in ax[1:]: a_.axvline(1,color='#888',ls=':',lw=.8); a_.set_xlim(0.4,50)
from matplotlib.lines import Line2D
ax[2].legend(handles=[Line2D([],[],color='k',lw=1,label='block-aligned'),Line2D([],[],color='k',lw=1,ls='--',label='through-block')],loc='upper right',fontsize=6)
ax[2].text(1.08,0.46,r'$T_{\rm on}$',fontsize=7,color='#555')
save(fig,'fig_edgeT')
