import pickle, glob, numpy as np, math
R={}
for f in glob.glob('dip/binder_*_s*.pkl'):
    L=int(f.split('_')[1]); R.setdefault(L,[]).append(pickle.load(open(f,'rb')))
Ls=sorted(R); TS=[r['T'] for r in R[Ls[0]][0]]
print('replicas:',{L:len(R[L]) for L in Ls})
def stat(L,i,key):
    v=np.array([run[i][key] for run in R[L]]); return v.mean(), v.std(ddof=1)/math.sqrt(len(v)) if len(v)>1 else 0
print('T      '+'  '.join(f'U(L={L})       ' for L in Ls)+'  2x (fit)      |<v>| L='+str(Ls[-1]))
rng=np.random.default_rng(0)
out=[]
for i,T in enumerate(TS):
    Us=[stat(L,i,'U') for L in Ls]
    # bootstrap over replicas for 2x
    xs=[]
    for b in range(400):
        m=[np.mean([R[L][k][i]['m2'] for k in rng.integers(0,len(R[L]),len(R[L]))]) for L in Ls]
        xs.append(-np.polyfit(np.log(Ls),np.log(m),1)[0])
    v=stat(Ls[-1],i,'vabs'); C=[stat(L,i,'C') for L in Ls]
    out.append(dict(T=T,U=Us,x2=(np.mean(xs),np.std(xs)),C=C,v=v))
    print(f'{("inf" if T>1e3 else f"{T:.2f}"):6s} '+'  '.join(f'{u:.3f}+-{e:.3f}' for u,e in Us)+f'  {np.mean(xs):.3f}+-{np.std(xs):.3f}  {v[0]:.3f}  C:'+' '.join(f'{c:.2f}' for c,e in C))
pickle.dump(dict(Ls=Ls,TS=TS,out=out,R=R),open('dip/binder_summary.pkl','wb'))
