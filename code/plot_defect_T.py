import pickle, numpy as np, math
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R=pickle.load(open('exact/pair_free_64.pkl','rb'))
fig,ax=plt.subplots(1,2,figsize=(14,5.2))
a=ax[0]
for k,c,lab,s in (('O','#6a1b9a','O pair (full vortices)',-2),('M','#d00000','monomer pair (half vortices)',-.5)):
    r=np.array([x for x,_ in R[k]]); z=np.array([y for _,y in R[k]])
    a.plot(r,z,'o',color=c,label=lab); rr=np.linspace(3,19,20); a.plot(rr,z[r==5][0]+s*np.log(rr/5),'--',color=c,lw=1,label=f'slope {s} (2x = {-s:g})')
a.set_xscale('log'); a.set_xlabel('separation r (cells)'); a.set_ylabel('ln [Z_pair(r)/Z]  (exact, 64×64 board)'); a.legend(fontsize=9); a.grid(alpha=.3,which='both')
a.set_title('(a) defect scaling dimensions: x_O = 1, x_monomer = 1/4\nx_O < 2 → free O defects are relevant (screen the order)',fontsize=11)
a=ax[1]; T=np.linspace(0.5,6,200)
for dO,c in ((4.4,'#08306b'),(2.0,'#4292c6'),(1.0,'#9ecae1')):
    a.semilogy(T,np.exp(dO/T),color=c,label=f'Δ_O = {dO} D per defect'+(' (dumbbell λ = 0.7)' if dO==4.4 else ''))
a.axvspan(0,0.45,color='#fde0dd'); a.text(0.47,2e3,'crystal\n(T < T_c ≈ 0.45 D)',fontsize=9)
a.set_xlim(0,6); a.set_ylim(1,1e4); a.set_xlabel('T (units D)'); a.set_ylabel('ξ / a₀ ~ 1/y = e^{Δ_O/T}'); a.legend(fontsize=9); a.grid(alpha=.3,which='both')
a.set_title('(b) screening length set by O fugacity (scaling, prefactor a₀ ~ 1 cell)\nmonomers stay bound by strings; below T_c the crystal confines O too',fontsize=11)
plt.tight_layout(); plt.savefig('outputs/defects_temperature.png',dpi=130); plt.savefig('outputs/defects_temperature.pdf'); print('ok')
