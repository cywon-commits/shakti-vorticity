import sys, time, pickle, random
from core import *
def chi(c): return 1 if black(c) else -1
def accumulate(samples_iter, cells):
    A={k:Counter() for k in ('short','conn','sig','tau')}; N=0; glob=[]
    for d in samples_iter:
        spin,free,ring,cf=analyse(d,cells); assert cf==0
        N+=1; rs=set(ring)
        for c in free: A['short'][c]+=1
        for c,s in spin.items():
            A['short'][c]+=1; A['conn'][c]+=1; A['sig'][c]+=s; A['tau'][c]+=s*chi(c)
        inner=[c for c in spin if c not in rs]
        glob.append(dict(nshort=len(spin)+len(free), nconn=len(spin), nring=len(ring), nfree=len(free),
                         M=sum(spin[c] for c in inner), Mtau=sum(spin[c]*chi(c) for c in inner)))
    return A,N,glob
def enum(n):
    cells,sh=aztec(n); cl=sorted(cells,key=lambda c:(c[1],c[0]))
    def rec(used,ds):
        free=next((c for c in cl if c not in used),None)
        if free is None: yield list(ds); return
        for dd in [(1,0),(0,1)]:
            o=(free[0]+dd[0],free[1]+dd[1])
            if o in cells and o not in used: yield from rec(used|{free,o},ds+[(free,o)])
    return rec(frozenset(),[]),cells
def sample(n,M,seed):
    cells,sh=aztec(n); rng=random.Random(seed)
    def it():
        for _ in range(M):
            yield [((a[0]+sh[0],a[1]+sh[1]),(b[0]+sh[0],b[1]+sh[1])) for a,b in shuffle(n,rng)]
    return it(),cells
if __name__=='__main__':
    n=int(sys.argv[1]); t=time.time()
    if sys.argv[2]=='exact': it,cells=enum(n)
    else: it,cells=sample(n,int(sys.argv[2]),n)
    A,N,glob=accumulate(it,cells)
    pickle.dump(dict(n=n,cells=cells,A=A,N=N,glob=glob,mode=sys.argv[2]),open(f'stats_{n}_{sys.argv[2]}.pkl','wb'))
    import statistics as st
    f=[g['nconn']/g['nshort'] for g in glob]
    print(n,sys.argv[2],'N',N,'time %.1fs'%(time.time()-t),'ring',glob[0]['nring'],
          'mean frac connected %.3f'%st.mean(f),'<M> %.3f'%st.mean(g['M'] for g in glob),'<Mtau> %.3f'%st.mean(g['Mtau'] for g in glob))
