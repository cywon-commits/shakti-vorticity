import pickle, glob, math, cmath, numpy as np
from analyse_nb import fourier, amp
out=[]
for f in sorted(glob.glob('nb/a_I_*.pkl')):
    D=pickle.load(open(f,'rb')); n=D['n']; sh=D['shift']
    res={}
    for Rw in (4,8):
        A=amp(fourier(D['vmean'],n,sh,Rw)); Ab=[amp(fourier(V,n,sh,Rw)) for V in D['vblocks']]
        res[Rw]=(A,np.std(Ab,ddof=1)/2)
    O=n+1; ec=D['emean'][O,O]; eb=[E[O,O] for E in D['eblocks']]
    out.append((n,res,ec,np.std(eb,ddof=1)/2))
out.sort()
for n,res,ec,ee in out: print('n=%3d  A(r<=4)=%.3f±%.3f  A(r<=8)=%.3f±%.3f  eps_c=%+.3f±%.3f'%(n,*res[4],*res[8],ec,ee))
def fit(xs,ys,es):
    x=np.log(xs); y=np.log(ys); w=np.array(ys)/np.array(es)
    p,c=np.polyfit(x,y,1,w=w,cov=True); return -p[0],math.sqrt(c[0,0])
for Rw in (4,8):
    R=[o for o in out if o[0]>=32]
    print('local window r<=%d: exponent %.3f ± %.3f (n=32-195)'%(Rw,*fit([o[0] for o in R],[o[1][Rw][0] for o in R],[o[1][Rw][1] for o in R])))
    for k in range(4):
        Rk=[o for o in R if o[0]%4==k]; print('    n≡%d: %.3f ± %.3f'%(k,*fit([o[0] for o in Rk],[o[1][Rw][0] for o in Rk],[o[1][Rw][1] for o in Rk])))
R=[o for o in out if o[0]%4 in (0,3) and o[0]>=32]
print('|eps(centre)| (n≡0,3): exponent %.3f ± %.3f'%fit([o[0] for o in R],[abs(o[2]) for o in R],[o[3] for o in R]))
pickle.dump(out,open('nb/local.pkl','wb'))
