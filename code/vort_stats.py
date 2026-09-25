import sys, time, pickle
from vort_tiling import *
from stats_run import enum, sample
def chi(c): return 1 if black(c) else -1
def run(it,cells):
    A={k:Counter() for k in ('v','vchi','vshort','nshort')}; N=0; glob=[]; frs=0; ringpat=set()
    for doms in it:
        I,adj,bad,own=signed_graph(doms); v,fr,comps=frustration(len(doms),adj)
        frs+= (fr>0 or bad>0 or comps!=1)
        spin,free,ring,_=analyse(doms,cells)
        sc={short_cell(*doms[k]):k for k in range(len(doms))}
        ref=v[sc[ring[0]]]
        ringpat.add(tuple(v[sc[c]]*ref for c in ring))
        M=0
        for k,d in enumerate(doms):
            s=v[k]*ref; M+=s
            for c in d: A['v'][c]+=s; A['vchi'][c]+=s*chi(c)
            scell=short_cell(*d); A['vshort'][scell]+=s; A['nshort'][scell]+=1
        glob.append(dict(M=M,ntet=len(doms)))
        N+=1
    return A,N,glob,frs,ringpat
if __name__=='__main__':
    n=int(sys.argv[1]); mode=sys.argv[2]; t=time.time()
    it,cells=enum(n) if mode=='exact' else sample(n,int(mode),n+2000)
    A,N,glob,frs,rp=run(it,cells)
    pickle.dump(dict(n=n,cells=cells,A=A,N=N,glob=glob,mode=mode),open(f'vstats_{n}_{mode}.pkl','wb'))
    import statistics as st
    Ms=[g['M']/g['ntet'] for g in glob]
    print(n,mode,'N',N,'bad',frs,'ring patterns',len(rp),'time %.0fs'%(time.time()-t),'<m> %.3f'%st.mean(Ms),'sd(m) %.3f'%st.pstdev(Ms))
