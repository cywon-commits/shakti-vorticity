import pickle, numpy as np, math, scipy.sparse as sp, scipy.sparse.linalg as spl
def green_square(L):
    # Dirichlet lattice Laplacian on interior vertices 1..L-1 of [0,L]^2 (boundary vertices at 0 and L)
    m=L-1; I=sp.identity(m); T=sp.diags([-1,2,-1],[-1,0,1],shape=(m,m))
    A=(sp.kron(T,I)+sp.kron(I,T)).tocsc()
    lu=spl.splu(A)
    return lambda p: lu.solve(np.eye(m*m)[:, (p[0]-1)*m+(p[1]-1)]).reshape(m,m), m
out={}
for L in (32,64,96):
    D=pickle.load(open(f'exact/corr_square_{L}.pkl','rb')); Gx=D['G']; p=(L//2,L//2)
    gf,m=green_square(L); Gp=gf(p)
    rs=list(range(1,D['Rw']+1)); pred=[]; ex=[]
    for r in rs:
        q=(p[0]+r,p[1]); Gq=gf(q)
        var=Gp[p[0]-1,p[1]-1]+Gq[q[0]-1,q[1]-1]-2*Gp[q[0]-1,q[1]-1]
        pred.append(math.exp(-(math.pi/2)*var)); ex.append(abs(Gx[(p,(r,0))]))
    out[L]=(rs,ex,pred)
    # amplitude from even r>=4
    ev=[i for i,r in enumerate(rs) if r%2==0 and r>=4]
    A=np.exp(np.mean([math.log(ex[i]/pred[i]) for i in ev]))
    print(L,'A=%.4f'%A,' ratio exact/(A*pred) at even r:',' '.join('%d:%.3f'%(rs[i],ex[i]/(A*pred[i])) for i in ev))
pickle.dump(out,open('exact/gff_pred.pkl','wb'))
