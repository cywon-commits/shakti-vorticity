import pickle
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
rows=pickle.load(open('hv_rows.pkl','rb'))
geoms=['Aztec I','Aztec II','Square, corner black','Square, corner white']
cols={0:'#1b6ca8',1:'#e0892b',2:'#3a9d5d',3:'#b0413e'}; mk={0:'o',1:'s',2:'^',3:'D'}
fig,ax=plt.subplots(1,4,figsize=(20,4.6),sharey=True)
for a,g in zip(ax,geoms):
    for k in range(4):
        R=[r for r in rows if r['geom']==g and r['size']%4==k]
        if R: a.plot([r['size'] for r in R],[r['asym'] for r in R],mk[k],color=cols[k],ms=6,label=f'size ≡ {k} (mod 4)')
    a.axhline(0,color='#bbb',lw=.8); a.set_title(g,fontsize=12); a.set_xlabel('n' if 'Aztec' in g else 'L'); a.grid(alpha=.25)
ax[0].set_ylabel('(A_H² − A_V²)/(A_H² + A_V²)\n+1: order only on horizontal plaquettes\n−1: only on vertical plaquettes')
ax[0].legend(frameon=False,fontsize=9); ax[0].set_ylim(-1.1,1.1)
plt.tight_layout(); plt.savefig('outputs/hv_split.png',dpi=180); plt.savefig('outputs/hv_split.pdf')
