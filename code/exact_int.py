"""Exact <eps(centre)> as a rational number via fraction-free (Bareiss) determinants with a real Kasteleyn sign."""
import sys, time
from fractions import Fraction
from exact_v import aztec_cells, Region, start_vertex, black
from vform2 import r4
def bareiss_det(M):
    M=[row[:] for row in M]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            for i in range(k+1,n):
                if M[i][k]!=0: M[k],M[i]=M[i],M[k]; sign=-sign; break
            else: return 0
        pk=M[k][k]; rowk=M[k]
        for i in range(k+1,n):
            mik=M[i][k]; rowi=M[i]
            if mik==0:
                for j in range(k+1,n): rowi[j]=rowi[j]*pk//prev
            else:
                for j in range(k+1,n): rowi[j]=(rowi[j]*pk-mik*rowk[j])//prev
            rowi[k]=0
        prev=pk
    return sign*M[n-1][n-1]
def real_K(R,flips):
    B,W=R.B,R.W; M=[[0]*len(W) for _ in B]
    for b in B:
        for di,dj in ((1,0),(-1,0),(0,1),(0,-1)):
            w=(b[0]+di,b[1]+dj)
            if w in R.wi:
                s=1 if dj==0 else (-1 if (min(b[1],w[1])+b[0])%2 else 1)   # vertical sign alternating with the column... checked below
                if (R.bi[b],R.wi[w]) in flips: s=-s
                M[R.bi[b]][R.wi[w]]=s
    return M
def flips_and_h0(R,start,hs,p):
    pts=R.path(start,p); h0=hs; flips=set()
    for a,b in zip(pts,pts[1:]):
        L,Rr=R.edge_cells(a,b); h0+= 1 if black(L) else -1
        if L in R.cells and Rr in R.cells:
            bc,wc=(L,Rr) if black(L) else (Rr,L); flips.add((R.bi[bc],R.wi[wc]))
    return flips,h0
for n in [int(x) for x in sys.argv[1:]]:
    t=time.time(); cells,sh=aztec_cells(n); R=Region(cells); s,hs=start_vertex(cells)
    Z=bareiss_det(real_K(R,set()))
    assert abs(Z)==2**(n*(n+1)//2), ('Kasteleyn sign wrong', Z)
    fl,h0=flips_and_h0(R,s,hs,(sh[0],0)); Zg=bareiss_det(real_K(R,fl))
    e0=1 if ((h0-r4((sh[0],0)))//4)%2==0 else -1
    val=Fraction(Zg,Z)*e0
    print(n, val, float(val), '(%.0fs)'%(time.time()-t), flush=True)
