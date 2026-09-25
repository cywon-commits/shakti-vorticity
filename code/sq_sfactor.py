"""MC structure factor of the vorticity field v on the Shakti plaquette centres (square L x L board)."""
import sys, time, pickle, math, numpy as np, io, contextlib
from square_mc import setup, initial_height, valid_edges, sweep, vfield
with contextlib.redirect_stdout(io.StringIO()):
    from core import pcen
L=int(sys.argv[1]); burn=int(sys.argv[2]); nm=int(sys.argv[3]); gap=int(sys.argv[4])
rng=np.random.default_rng(L)
X,Y,sgn,R=setup(L,0,0); h=initial_height(L,0,0,X,Y,sgn,R); assert valid_edges(h,sgn)
P=np.array([[pcen((i,j)) for j in range(L)] for i in range(L)])      # L x L x 2
gx=np.rint(2*P[...,0]).astype(int); gy=np.rint(2*P[...,1]).astype(int); gx-=gx.min(); gy-=gy.min()
M=1<<int(math.ceil(math.log2(max(gx.max(),gy.max())+1)))
t=time.time()
for _ in range(burn): sweep(h,sgn,rng,L)
acc=np.zeros((M,M)); vsum=np.zeros((L,L))
for m in range(nm):
    for _ in range(gap): sweep(h,sgn,rng,L)
    v,_,_=vfield(h,X,Y,sgn,L); g=np.zeros((M,M)); g[gx,gy]=v
    acc+=np.abs(np.fft.fft2(g))**2; vsum+=v
g=np.zeros((M,M)); g[gx,gy]=vsum/nm; Sb=np.abs(np.fft.fft2(g))**2/(L*L)
S=acc/nm/(L*L); q=np.fft.fftfreq(M,d=0.5)*2*math.pi
pickle.dump(dict(L=L,M=M,q=q,S=S,Sbragg=Sb,N=nm,burn=burn,gap=gap),open(f'sq/sf_{L}.pkl','wb'))
print(L,'M',M,'done %.0fs'%(time.time()-t))
