"""Dipolar Metropolis (ground-state manifold, point dipoles): independent annealing replicas for error bars,
Binder cumulant of the (pi/2,pi/2) crystal order parameter, and the effective exponent 2x from <|Psi|^2>.
Usage: python dip_binder.py L seed neq nsamp gap"""
import sys, math, pickle, time, numpy as np
from dip_mc import *
from square_mc import vfield
L=int(sys.argv[1]); seed=int(sys.argv[2]); neq,nsamp,gap=map(int,sys.argv[3:6])
TS=[1e9,1.0,0.8,0.7,0.6,0.55,0.5,0.47,0.45,0.43,0.4,0.35]
D=pickle.load(open(f'dip/setup_{L}_point.pkl','rb'))
X,Y,sgn,R=setup(L,0,0); h=initial_height(L,0,0,X,Y,sgn,R).astype(np.int64)
px,mx,py,my=[np.ascontiguousarray(sgn[k].astype(np.int64)) for k in ('px','mx','py','my')]
s=np.zeros(len(D['E']),np.int8); init_spins(h,px,py,L,D['Ti'],D['Ts'],D['Tl'],D['sd'],s)
Et=np.array([total_energy(s,D['NBi'],D['NBJ'],D['NBl'])])
lo,hi=L//4,3*L//4; x=np.arange(L)[:,None]; y=np.arange(L)[None,:]; ph=np.exp(1j*math.pi/2*(x+y))
dummy=np.zeros(2); res=[]; sd=seed*1000
for T in TS:
    t=time.time()
    run(h,px,mx,py,my,L,D['Ti'],D['Ts'],D['Tl'],D['sd'],s,D['NBi'],D['NBJ'],D['NBl'],T,neq,sd,False,1,dummy,dummy,dummy,Et); sd+=1
    m2=[];c2=[];c4=[];E=[];va=[]
    for k in range(nsamp):
        run(h,px,mx,py,my,L,D['Ti'],D['Ts'],D['Tl'],D['sd'],s,D['NBi'],D['NBJ'],D['NBl'],T,gap,sd,False,1,dummy,dummy,dummy,Et); sd+=1
        w=np.exp(1j*math.pi*h[lo:hi,lo:hi]/4); m2.append(abs(w.mean())**2)
        v,Hd,Vd=vfield(h,X,Y,sgn,L); c=(ph*(Hd.astype(float)-Vd.astype(float))).mean()
        c2.append(abs(c)**2); c4.append(abs(c)**4); E.append(Et[0]); va.append(v[lo:hi,lo:hi].astype(float))
    E=np.array(E)
    r=dict(T=T,m2=np.mean(m2),c2=np.mean(c2),c4=np.mean(c4),U=1-np.mean(c4)/(2*np.mean(c2)**2),E=E.mean()/L**2,C=E.var()/T**2/L**2,vabs=np.mean(np.abs(np.mean(va,axis=0))))
    res.append(r); print(f'L={L} seed={seed} T={T:.3f} m2={r["m2"]:.4f} U={r["U"]:.3f} C={r["C"]:.3f} |<v>|={r["vabs"]:.3f} {time.time()-t:.0f}s',flush=True)
pickle.dump(res,open(f'dip/binder_{L}_s{seed}.pkl','wb'))
