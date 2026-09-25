import pickle, math, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
fig,ax=plt.subplots(2,2,figsize=(12,11))
for row,n in enumerate((20,40)):
    D=pickle.load(open(f'sfactor_{n}.pkl','rb')); q=D['q']; S=D['S']; Sb=D['Sbragg']
    ext=(q[0]/math.pi,q[-1]/math.pi,q[0]/math.pi,q[-1]/math.pi)
    for col,(Z,lab) in enumerate(((S,'total S(q)'),(Sb,'static part |⟨m⟩_q|² (Bragg)'))):
        a=ax[row][col]; Zp=np.maximum(Z,1e-4)
        im=a.imshow(Zp.T,origin='lower',extent=ext,norm=LogNorm(vmin=max(Zp.max()*1e-4,1e-4),vmax=Zp.max()),cmap='magma')
        a.set_xlim(-2,2); a.set_ylim(-2,2); a.set_xlabel('q_x / π'); a.set_ylabel('q_y / π'); a.set_title(f'n = {n}: {lab}',fontsize=11)
        fig.colorbar(im,ax=a,shrink=.8)
    # strongest static peaks away from q=0
    Q=np.array(np.meshgrid(q,q,indexing='ij')); mask=(np.hypot(Q[0],Q[1])>0.3)
    idx=np.argsort((Sb*mask).ravel())[::-1][:12]
    pk=sorted({(round(q[i//len(q)]/math.pi,2),round(q[i%len(q)]/math.pi,2)) for i in idx})
    print(n,'strongest static peaks (q/π):',pk)
fig.suptitle('Predicted magnetic structure factor of the island moments (Aztec diamond, placement I, transverse projection)',fontsize=12)
plt.tight_layout(rect=(0,0,1,.96)); plt.savefig('outputs/structure_factor.png',dpi=150); plt.savefig('outputs/structure_factor.pdf')
