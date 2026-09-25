"""Metropolis on domino tilings (height flips) with the dipolar energy of the tetromino-template spin configuration."""
import sys, math, pickle, time, numpy as np, numba as nb
from square_mc import setup, initial_height, valid_edges
@nb.njit(cache=True)
def epsv(hv): return 1 if (hv//4)%2==0 else -1
@nb.njit(cache=True)
def local_field(k,s,NBi,NBJ,NBl):
    f=0.0
    for m in range(NBl[k]): f+=NBJ[k,m]*s[NBi[k,m]]
    return f
@nb.njit(cache=True)
def total_energy(s,NBi,NBJ,NBl):
    E=0.0
    for k in range(s.shape[0]): E+=0.5*s[k]*local_field(k,s,NBi,NBJ,NBl)
    return E
@nb.njit(cache=True)
def init_spins(h,px,py,L,Ti,Ts,Tl,sd,s):
    for i in range(L):
        for j in range(L):
            if i<L-1 and (h[i+1,j+1]-h[i+1,j])*py[i+1,j]==-3:
                d=i*L+j; v=sd[d]*epsv(h[i,j])
                for m in range(Tl[d]): s[Ti[d,m]]=v*Ts[d,m]
            if j<L-1 and (h[i+1,j+1]-h[i,j+1])*px[i,j+1]==-3:
                d=L*L+i*L+j; v=sd[d]*epsv(h[i,j])
                for m in range(Tl[d]): s[Ti[d,m]]=v*Ts[d,m]
@nb.njit(cache=True)
def run(h,px,mx,py,my,L,Ti,Ts,Tl,sd,s,NBi,NBJ,NBl,T,nsweep,seed,meas,rmax,acc_dh,acc_dh2,acc_n,Etrack):
    np.random.seed(seed)
    ni=s.shape[0]; mark=np.zeros(ni,np.int64); snew=np.zeros(ni,np.int8)
    C=np.zeros(128,np.int64); nacc=0; E=Etrack[0]
    nin=(L-1)*(L-1)
    for sw in range(nsweep):
        for it in range(nin):
            a=1+np.random.randint(L-1); b=1+np.random.randint(L-1)
            prop=h[a,b]+(4 if np.random.random()<0.5 else -4)
            ok=True
            d1=(h[a+1,b]-prop)*px[a,b]; d2=(h[a-1,b]-prop)*mx[a,b]; d3=(h[a,b+1]-prop)*py[a,b]; d4=(h[a,b-1]-prop)*my[a,b]
            if not ((d1==1 or d1==-3) and (d2==1 or d2==-3) and (d3==1 or d3==-3) and (d4==1 or d4==-3)): continue
            upcross=(h[a,b+1]-h[a,b])*py[a,b]==-3
            if upcross:   # H pair -> V pair
                n1=L*L+(a-1)*L+(b-1); n2=L*L+a*L+(b-1)
            else:
                n1=(a-1)*L+(b-1); n2=(a-1)*L+b
            nc=0
            for d in (n1,n2):
                pi=(d%(L*L))//L; pj=d%L
                v=sd[d]*epsv(h[pi,pj])
                for m in range(Tl[d]):
                    k=Ti[d,m]; val=v*Ts[d,m]
                    if val!=s[k] and mark[k]==0:
                        mark[k]=1; snew[k]=val; C[nc]=k; nc+=1
            dE=0.0
            for x in range(nc):
                k=C[x]; dk=snew[k]-s[k]
                dE+=dk*local_field(k,s,NBi,NBJ,NBl)
            for x in range(nc):
                k=C[x]; dk=snew[k]-s[k]
                for m in range(NBl[k]):
                    l=NBi[k,m]
                    if mark[l]==1 and l>k: dE+=NBJ[k,m]*dk*(snew[l]-s[l])
            if dE<=0 or np.random.random()<math.exp(-dE/T):
                h[a,b]=prop; E+=dE; nacc+=1
                for x in range(nc): s[C[x]]=snew[C[x]]
            for x in range(nc): mark[C[x]]=0
        if meas:
            lo=L//4; hi=3*L//4
            for r in range(1,rmax+1):
                for x in range(lo,hi-r+1):
                    for y in range(lo,hi+1):
                        dh=h[x+r,y]-h[x,y]; acc_dh[r]+=dh; acc_dh2[r]+=dh*dh
                        dh=h[y,x+r]-h[y,x]; acc_dh[r]+=dh; acc_dh2[r]+=dh*dh
                        acc_n[r]+=2
    Etrack[0]=E
    return nacc
def simulate(L,model,Ts_list,neq,nmeas,seed=1,tag=''):
    D=pickle.load(open(f'dip/setup_{L}_{model}.pkl','rb'))
    X,Y,sgn,R=setup(L,0,0); h=initial_height(L,0,0,X,Y,sgn,R).astype(np.int64); assert valid_edges(h,sgn)
    px,mx,py,my=[np.ascontiguousarray(sgn[k].astype(np.int64)) for k in ('px','mx','py','my')]
    s=np.zeros(len(D['E']),np.int8); init_spins(h,px,py,L,D['Ti'],D['Ts'],D['Tl'],D['sd'],s)
    Et=np.array([total_energy(s,D['NBi'],D['NBJ'],D['NBl'])])
    rmax=L//4; out=[]
    for T in Ts_list:
        a=np.zeros(rmax+1); a2=np.zeros(rmax+1); an=np.zeros(rmax+1)
        t=time.time()
        run(h,px,mx,py,my,L,D['Ti'],D['Ts'],D['Tl'],D['sd'],s,D['NBi'],D['NBJ'],D['NBl'],T,neq,seed,False,rmax,a,a2,an,Et)
        nacc=run(h,px,mx,py,my,L,D['Ti'],D['Ts'],D['Tl'],D['sd'],s,D['NBi'],D['NBJ'],D['NBl'],T,nmeas,seed+7,True,rmax,a,a2,an,Et)
        Echeck=total_energy(s,D['NBi'],D['NBJ'],D['NBl'])
        var=a2[1:]/an[1:]-(a[1:]/an[1:])**2
        rs=np.arange(1,rmax+1); sel=(rs>=3)
        kappa=np.polyfit(np.log(rs[sel]),var[sel],1)[0]
        out.append(dict(T=T,var=var,kappa=kappa,x=math.pi**2*kappa/64,acc=nacc/(nmeas*(L-1)**2),E=Et[0]/L**2,Echeck=Echeck/L**2,h=h.copy(),s=s.copy()))
        print(f'{tag} L={L} {model} T={T:.3f}: kappa={kappa:.3f} (free 1.621)  x={math.pi**2*kappa/64:.3f}  acc={nacc/(nmeas*(L-1)**2):.3f}  E/L^2={Et[0]/L**2:.4f} (check {Echeck/L**2:.4f})  {time.time()-t:.0f}s',flush=True)
        seed+=13
    return out
if __name__=='__main__':
    L=int(sys.argv[1]); model=sys.argv[2]; neq=int(sys.argv[3]); nmeas=int(sys.argv[4])
    Ts_list=[float(x) for x in sys.argv[5].split(',')]
    out=simulate(L,model,Ts_list,neq,nmeas,tag=sys.argv[6] if len(sys.argv)>6 else '')
    pickle.dump(out,open(f'dip/run_{L}_{model}_{sys.argv[6] if len(sys.argv)>6 else "x"}.pkl','wb'))
