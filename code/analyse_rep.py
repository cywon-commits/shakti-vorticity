import pickle, glob, math, numpy as np
from analyse_square import amp
res={}
for cls,sh in (('black','00'),('white','10')):
    for L in (16,24,32,40,48):
        A=[amp(pickle.load(open(f,'rb')))[0] for f in sorted(glob.glob(f'sqrep/sq_{L}_{sh}_r*.pkl'))]
        full=[r for r in pickle.load(open('sq/rows.pkl','rb')) if r['cls']==cls and r['L']==L][0]['amp']
        m=np.mean(A); e=np.std(A,ddof=1)/math.sqrt(len(A))
        res[(cls,L)]=(m,e,full,len(A))
        print('%s L=%d  replicas %d: A=%.3f ± %.3f   (long run: %.3f)'%(cls,L,len(A),m,e,full))
for cls in ('black','white'):
    Ls=[L for L in (16,24,32,40,48)]; y=np.log([res[(cls,L)][0] for L in Ls]); w=[res[(cls,L)][0]/res[(cls,L)][1] for L in Ls]
    p,c=np.polyfit(np.log(Ls),y,1,w=w,cov=True); print(cls,'exponent %.3f ± %.3f (L=16-48, replica errors)'%(-p[0],math.sqrt(c[0,0])))
pickle.dump(res,open('sqrep/res.pkl','wb'))
