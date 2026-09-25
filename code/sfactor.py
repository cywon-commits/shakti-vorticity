"""Predicted island-moment structure factor (Aztec diamond, placement I), via FFT on a half-lattice grid."""
import sys, math, pickle, numpy as np, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    from vort import tet_of, template_for
from aztec_nb import shuffle, height, SH, SV, R4
n=int(sys.argv[1]); NS=int(sys.argv[2]); shx,shy=(n+1)%2,0; O=n+1
cache={}
def tmpl(d):
    if d not in cache:
        rs_,rl,s=tet_of(d); tp,fr=template_for(rs_,rl); cache[d]={e:sg for e,(sg,k) in tp.items() if k!='ext'}
    return cache[d]
G=None; acc=None; mx_sum=None; my_sum=None
for t in range(NS):
    T=shuffle(n,50000+t); h=height(T,n,shx,shy,R4); spins={}
    for i in range(-n,n):
        for j in range(-n,n):
            tp=T[i+O,j+O]
            if tp==0: continue
            x,y=i+shx,j+shy; hv=h[i+O,j+O]; e=1 if ((hv-(hv%4))//4)%2==0 else -1
            if tp in (1,2): d=((x,y),(x+1,y)); v=SH[x%4,y%4]*e
            else: d=((x,y),(x,y+1)); v=SV[x%4,y%4]*e
            for edge,sg in tmpl(d).items(): spins[edge]=v*sg
    if G is None:
        E=sorted(spins); xs=[a[0]+b[0] for a,b in E]; ys=[a[1]+b[1] for a,b in E]
        x0,y0=min(xs),min(ys); M=256
        G=(np.array(xs)-x0, np.array(ys)-y0, E)
        acc=np.zeros((M,M)); mx_sum=np.zeros((M,M)); my_sum=np.zeros((M,M))
    gx,gy,E=G; mx=np.zeros((M,M)); my=np.zeros((M,M))
    for k,(a,b) in enumerate(E):
        s=spins[(a,b)]; mx[gx[k],gy[k]]=(b[0]-a[0])*s; my[gx[k],gy[k]]=(b[1]-a[1])*s
    Fx=np.fft.fftshift(np.fft.fft2(mx)); Fy=np.fft.fftshift(np.fft.fft2(my))
    q=np.fft.fftshift(np.fft.fftfreq(M,d=0.5))*2*math.pi; QX,QY=np.meshgrid(q,q,indexing='ij'); q2=QX**2+QY**2; q2[q2==0]=1
    acc+=(np.abs(Fx)**2+np.abs(Fy)**2-np.abs(QX*Fx+QY*Fy)**2/q2)/len(E)
    mx_sum+=mx; my_sum+=my
mxa=mx_sum/NS; mya=my_sum/NS
Fx=np.fft.fftshift(np.fft.fft2(mxa)); Fy=np.fft.fftshift(np.fft.fft2(mya))
Sb=(np.abs(Fx)**2+np.abs(Fy)**2-np.abs(QX*Fx+QY*Fy)**2/q2)/len(E)
pickle.dump(dict(n=n,NS=NS,q=q,S=acc/NS,Sbragg=Sb,grid=(gx,gy),mx=mxa,my=mya),open(f'sfactor_{n}.pkl','wb'))
print('done',n,NS,'islands',len(E))
