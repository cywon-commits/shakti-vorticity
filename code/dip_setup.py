"""Precompute arrays for dipolar Metropolis on domino tilings of an L x L chess board (black corner at (0,0))."""
import math, numpy as np, io, contextlib, pickle, sys
with contextlib.redirect_stdout(io.StringIO()):
    from vort import template_for, key
    from core import rect, short_cell
from vform2 import S_TAB, r4
def build(L,Rc=6.0,model='point',lam=0.7):
    doms=[]
    for i in range(L):
        for j in range(L):
            doms.append(((i,j),(i+1,j)) if i<L-1 else None)
    for i in range(L):
        for j in range(L):
            doms.append(((i,j),(i,j+1)) if j<L-1 else None)
    isl={}; T=[]; sd=np.zeros(len(doms),np.int8)
    for k,d in enumerate(doms):
        if d is None: T.append([]); continue
        s=short_cell(*d); l=d[1] if s==d[0] else d[0]
        tp,fr=template_for(list(rect(s)),list(rect(l)))
        row=[]
        for e,(sg,kind) in tp.items():
            if e not in isl: isl[e]=len(isl)
            row.append((isl[e],sg))
        T.append(row)
        typ='H' if d[1][0]!=d[0][0] else 'V'; p=d[0]; sd[k]=S_TAB[(typ,p[0]%4,p[1]%4)]
    M=max(len(r) for r in T); nd=len(doms)
    Ti=-np.ones((nd,M),np.int64); Ts=np.zeros((nd,M),np.int8); Tl=np.zeros(nd,np.int64)
    for k,r in enumerate(T):
        Tl[k]=len(r)
        for m,(a,b) in enumerate(r): Ti[k,m]=a; Ts[k,m]=b
    E=sorted(isl,key=lambda e:isl[e]); ni=len(E)
    pos=np.array([[(a[0]+b[0])/2,(a[1]+b[1])/2] for a,b in E]); dirv=np.array([[b[0]-a[0],b[1]-a[1]] for a,b in E],float)
    # couplings within cutoff
    nb=[[] for _ in range(ni)]
    for k in range(ni):
        d=pos-pos[k]; r=np.hypot(d[:,0],d[:,1]); idx=np.where((r>0)&(r<=Rc))[0]
        for l in idx:
            if model=='point':
                rv=d[l]/r[l]; J=(dirv[k]@dirv[l]-3*(dirv[k]@rv)*(dirv[l]@rv))/r[l]**3
            else:   # dumbbell: charges +-1/lam at the ends of a needle of length lam, energy sum q q'/r
                J=0.0
                for sa in (1,-1):
                    for sb in (1,-1):
                        pa=pos[k]+sa*lam/2*dirv[k]; pb=pos[l]+sb*lam/2*dirv[l]; J+=sa*sb/np.hypot(*(pa-pb))
                J/=lam*lam
            nb[k].append((l,J))
    K=max(len(x) for x in nb); NBi=-np.ones((ni,K),np.int64); NBJ=np.zeros((ni,K)); NBl=np.zeros(ni,np.int64)
    for k,x in enumerate(nb):
        NBl[k]=len(x)
        for m,(l,J) in enumerate(x): NBi[k,m]=l; NBJ[k,m]=J
    return dict(L=L,Ti=Ti,Ts=Ts,Tl=Tl,sd=sd,NBi=NBi,NBJ=NBJ,NBl=NBl,pos=pos,dirv=dirv,E=E,model=model,Rc=Rc)
if __name__=='__main__':
    L=int(sys.argv[1]); model=sys.argv[2] if len(sys.argv)>2 else 'point'
    D=build(L,model=model); pickle.dump(D,open(f'dip/setup_{L}_{model}.pkl','wb')); print(L,model,'islands',len(D['E']),'max nb',D['NBi'].shape[1])
