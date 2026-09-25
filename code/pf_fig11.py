from pubstyle import *
import pickle, random, math, numpy as np
src=open('plot_mono.py').read().split("L=14; board")[0]
src=src.replace("import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt","")
src=src.replace("ax.plot([C[i][0],C[j][0]],[C[i][1],C[j][1]],color='#d00000',lw=3.5,zorder=9,solid_capstyle='round')","ax.plot([C[i][0],C[j][0]],[C[i][1],C[j][1]],color='#d00000',lw=1.8,zorder=9,solid_capstyle='round')")
src=src.replace("outline(ax,s+l); C.append(cen(s+l))","outline(ax,s+l,lw=.5); C.append(cen(s+l))")
src=src.replace("outline(ax,sq,lw=1.6); c=cen(sq); ax.plot(*c,'k*',ms=11,zorder=8)","outline(ax,sq,lw=.8); c=cen(sq); ax.plot(*c,'k*',ms=5,zorder=8)")
src=src.replace("ax.text(c[0],c[1],'O',ha='center',va='center',color='white',fontsize=8,fontweight='bold')","ax.text(c[0],c[1],'O',ha='center',va='center',color='white',fontsize=5,fontweight='bold')")
src=src.replace("    ax.set_title(f'{title}\\nminimum # violated couplings = {f}',fontsize=11)","    ax.set_title(title+f' [{f}]',loc='left',pad=2)")
exec(src)
fig=plt.figure(figsize=(TWO,4.0))
gs=fig.add_gridspec(2,3,height_ratios=[1.2,1],wspace=0.32,hspace=0.28,left=0.06,right=0.99,top=0.95,bottom=0.09)
L=14; board={(i,j) for i in range(L) for j in range(L)}; rng=random.Random(11)
a,b=(3,7),(10,7)
ax1=fig.add_subplot(gs[0,0]); tiles=random_tiling(board-{a,b},rng); draw(ax1,tiles,[a,b],'(a) monomer pair')
Oa=tuple(sorted((a,partner(a)))); bb=(9,8); Ob=tuple(sorted((bb,partner(bb))))
ax2=fig.add_subplot(gs[0,1]); tiles=random_tiling(board-set(Oa)-set(Ob),rng)+[('O',)+Oa,('O',)+Ob]; draw(ax2,tiles,[],'(b) diagonal-tile pair')
ax3=fig.add_subplot(gs[0,2]); tiles=random_tiling(board-{(6,6),(0,13)},rng); draw(ax3,tiles,[(6,6),(0,13)],'(c) string to the edge')
R=pickle.load(open('mono/tile_scan.pkl','rb'))
a4=fig.add_subplot(gs[1,0]); ds=sorted({d for d,_ in R['mono']})
for d,fs in R['mono']: a4.plot([d]*len(fs),fs,'o',color='#d00000',alpha=.25,ms=2.5,mew=0)
a4.plot(ds,[np.mean(sum([fs for d2,fs in R['mono'] if d2==d],[])) for d in ds],'-o',color='#d00000',ms=3,lw=.9,label='$\\Delta h=\\pm4$ (monomers)')
a4.plot([dx+1 for dx,fs in R['O']],[np.mean(fs) for dx,fs in R['O']],'s',color='k',ms=3.5,label='$\\Delta h=\\pm8$ (diagonal tiles)')
a4.set_xlabel('separation (cells)'); a4.set_ylabel('violated couplings'); a4.legend(loc='upper left',fontsize=6.5); a4.set_title('(d)',loc='left',pad=2)
a5=fig.add_subplot(gs[1,1])
for k,D in enumerate((2,3,5)):
    res=pickle.load(open(f'mono/path24_D{D}.pkl','rb')); keep=[];rm=0
    for r,e,x in res:
        if r>rm+1e-9 and r<13.3: keep.append((r,e)); rm=r
        else: break
    a5.plot([k_[0] for k_ in keep],[k_[1] for k_ in keep],'-',marker=MK[k],color=CAT[k],ms=3,lw=.9,label=f'$\\Delta={D}$')
a5.set_xlabel('separation $r$ ($a$)'); a5.set_ylabel('excess energy ($J$)'); a5.legend(loc='upper left'); a5.set_title('(e)',loc='left',pad=2)
P=pickle.load(open('exact/pair_free_64.pkl','rb'))
a6=fig.add_subplot(gs[1,2])
for k,c,lab,sl,m in (('O','k','$\\Delta h=\\pm8$',-2,'s'),('M','#d00000','$\\Delta h=\\pm4$',-.5,'o')):
    r=np.array([x for x,_ in P[k]]); z=np.array([y for _,y in P[k]])
    a6.plot(r,z,m,color=c,ms=3,label=lab); rr=np.linspace(3,19,20); a6.plot(rr,z[r==5][0]+sl*np.log(rr/5),'--',color=c,lw=.7)
a6.set_xscale('log'); a6.set_xlabel('separation $r$ (cells)'); a6.set_ylabel('$\\ln[Z_{\\rm pair}(r)/Z]$'); a6.legend(loc='lower left')
a6.text(6,-2.0,'slope $-1/2$',color='#d00000',fontsize=6.5); a6.text(9,-6.3,'slope $-2$',fontsize=6.5); a6.set_title('(f)',loc='left',pad=2)
import matplotlib.ticker as mt; a6.set_xticks([1,3,10]); a6.set_xticklabels(['1','3','10']); a6.xaxis.set_minor_formatter(mt.NullFormatter())
save(fig,'fig11_strings')
