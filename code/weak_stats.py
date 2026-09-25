import sys, time, pickle, random
sys_argv=sys.argv[:]; sys.argv=['x','1']
exec(open('weak_variants.py').read().split("n=int(sys.argv[1])")[0])
sys.argv=sys_argv
from stats_run import enum, sample
def chi(c): return 1 if black(c) else -1
def acc(it,cells,variant):
    A={k:Counter() for k in ('short','conn','sig','tau','ferro')}; N=0; glob=[]; frs=0
    for d in it:
        S,adj,ring,ferro=build_v(d,cells,variant)
        spin,comps,f=solve(S,adj,ring); frs+=bool(f); N+=1
        rs=set(ring)
        for s in S: A['short'][s]+=1
        for e in ferro:
            for x in e: A['ferro'][x]+=1
        for c,s in spin.items(): A['conn'][c]+=1; A['sig'][c]+=s; A['tau'][c]+=s*chi(c)
        inner=[c for c in spin if c not in rs]
        glob.append(dict(nshort=len(S),nconn=len(spin),nferro=len(ferro),Mtau=sum(spin[c]*chi(c) for c in inner),
                         ncomp=len(comps),maxfree=max([len(c) for c in comps],default=0)))
    return A,N,glob,frs
if __name__=='__main__':
    n=int(sys.argv[1]); mode=sys.argv[2]; variant=sys.argv[3] if len(sys.argv)>3 else 'one_isolated'; t=time.time()
    it,cells=enum(n) if mode=='exact' else sample(n,int(mode),n+1000)
    A,N,glob,frs=acc(it,cells,variant)
    pickle.dump(dict(n=n,cells=cells,A=A,N=N,glob=glob,mode=mode),open(f'wstats_{variant}_{n}_{mode}.pkl','wb'))
    import statistics as st
    print(n,mode,variant,'N',N,'frustrated',frs,'time %.0fs'%(time.time()-t),
      'frac fixed %.3f'%st.mean(g['nconn']/g['nshort'] for g in glob),'<Mtau> %.2f'%st.mean(g['Mtau'] for g in glob),
      '<#free comps> %.2f'%st.mean(g['ncomp'] for g in glob),'<#ferro> %.1f'%st.mean(g['nferro'] for g in glob))
