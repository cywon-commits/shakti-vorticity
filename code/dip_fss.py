import sys, math, pickle, time, numpy as np
from dip_mc import *
from square_mc import vfield
def fss(L,model,Ts_list,neq,nsamp,gap,seed=3):
    D=pickle.load(open(f'dip/setup_{L}_{model}.pkl','rb'))
    X,Y,sgn,R=setup(L,0,0); h=initial_height(L,0,0,X,Y,sgn,R).astype(np.int64)
    px,mx,py,my=[np.ascontiguousarray(sgn[k].astype(np.int64)) for k in ('px','mx','py','my')]
    s=np.zeros(len(D['E']),np.int8); init_spins(h,px,py,L,D['Ti'],D['Ts'],D['Tl'],D['sd'],s)
    Et=np.array([total_energy(s,D['NBi'],D['NBJ'],D['NBl'])])
    lo,hi=L//4,3*L//4; x=np.arange(L)[:,None]; y=np.arange(L)[None,:]
    ph=np.exp(1j*math.pi/2*(x+y))           # (pi/2,pi/2) crystal wave
    dummy=np.zeros(2); res=[]
    for T in Ts_list:
        run(h,px,mx,py,my,L,D['Ti'],D['Ts'],D['Tl'],D['sd'],s,D['NBi'],D['NBJ'],D['NBl'],T,neq,seed,False,1,dummy,dummy,dummy,Et); seed+=1
        m2=[];c2=[];E=[];vm=[]
        for k in range(nsamp):
            run(h,px,mx,py,my,L,D['Ti'],D['Ts'],D['Tl'],D['sd'],s,D['NBi'],D['NBJ'],D['NBl'],T,gap,seed,False,1,dummy,dummy,dummy,Et); seed+=1
            w=np.exp(1j*math.pi*h[lo:hi,lo:hi]/4); m2.append(abs(w.mean())**2)
            v,Hd,Vd=vfield(h,X,Y,sgn,L); c=(ph*(Hd.astype(float)-Vd.astype(float)))[lo:hi,lo:hi].mean(); c2.append(abs(c)**2)
            E.append(Et[0]); vm.append(v[lo:hi,lo:hi].copy())
        E=np.array(E); r=dict(T=T,m2=np.mean(m2),m2err=np.std(m2)/math.sqrt(nsamp),c2=np.mean(c2),E=E.mean()/L**2,C=E.var()/T**2/L**2,vmap=np.mean(vm,axis=0))
        res.append(r); print(f'L={L} T={T:.3f} <|Psi|^2>={r["m2"]:.4f}  crystal^2={r["c2"]:.4f}  E/L2={r["E"]:.4f}  C={r["C"]:.3f}',flush=True)
    return res
if __name__=='__main__':
    L=int(sys.argv[1]); model=sys.argv[2]
    Ts=[float(x) for x in sys.argv[3].split(',')]; neq,nsamp,gap=int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6])
    res=fss(L,model,Ts,neq,nsamp,gap); pickle.dump(res,open(f'dip/fss_{L}_{model}.pkl','wb'))
