import pickle, sys, numpy as np, math, cmath
from square_mc import setup, vfield
def diag(h,L,X,Y,sgn):
    v,Hd,Vd=vfield(h,X,Y,sgn,L)   # Hd[x,y]: horizontal domino with left cell (x,y)
    lo,hi=L//4,3*L//4; W=slice(lo,hi)
    x=np.arange(L)[:,None]; y=np.arange(L)[None,:]
    nH=Hd.astype(float); nV=Vd.astype(float)
    colx=np.mean(((-1)**x*nH)[W,W])*4; coly=np.mean(((-1)**y*nV)[W,W])*4
    st=np.mean(((-1)**(x+y)*(nH-nV))[W,W])*2
    fracH=nH[W,W].sum()/(nH[W,W].sum()+nV[W,W].sum())
    # flippable plaquettes
    fl=0; tot=0
    for a in range(lo,hi):
        for b in range(lo,hi):
            tot+=1
            if (Hd[a,b] and Hd[a,b+1]) or (Vd[a,b] and Vd[a+1,b]): fl+=1
    tilt=((h[hi,lo:hi]-h[lo,lo:hi]).mean()/(hi-lo), (h[lo:hi,hi]-h[lo:hi,lo]).mean()/(hi-lo))
    psi=np.mean(np.exp(1j*math.pi*h[lo:hi,lo:hi]/4))
    vm=np.mean(v[W,W])
    return dict(colx=colx,coly=coly,stag=st,fracH=fracH,flip=fl/tot,tilt=tilt,psi=abs(psi),vmean=vm)
if __name__=='__main__':
    out=pickle.load(open(sys.argv[1],'rb')); L=int(sys.argv[2]); X,Y,sgn,R=setup(L,0,0)
    for o in out:
        d=diag(o['h'],L,X,Y,sgn)
        print('T=%-8g col(%.2f,%.2f) stag %.2f fracH %.2f flippable %.2f tilt(%.2f,%.2f) |mean e^{iπh/4}| %.2f'%(o['T'],d['colx'],d['coly'],d['stag'],d['fracH'],d['flip'],*d['tilt'],d['psi']))
