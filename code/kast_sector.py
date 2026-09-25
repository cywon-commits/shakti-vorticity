"""Sector (winding) distribution of domino tilings of an L x L board with a central H x H hole, via Kasteleyn + flux."""
import sys, numpy as np, math
def annulus(L,H):
    o=(L-H)//2; cell=np.ones((L,L),bool); cell[o:o+H,o:o+H]=False; return cell,o
def kast(L,H,theta):
    cell,o=annulus(L,H)
    B=[(i,j) for i in range(L) for j in range(L) if cell[i,j] and (i+j)%2==0]
    W=[(i,j) for i in range(L) for j in range(L) if cell[i,j] and (i+j)%2==1]
    bi={c:k for k,c in enumerate(B)}; wi={c:k for k,c in enumerate(W)}
    K=np.zeros((len(B),len(W)),complex)
    jc=L//2
    for (i,j) in B:
        for di,dj in ((1,0),(-1,0),(0,1),(0,-1)):
            q=(i+di,j+dj)
            if q in wi:
                wgt=1.0 if dj==0 else 1j
                # flux on vertical bonds crossing the line y = jc to the right of the hole
                if dj!=0 and i>=o+H and min(j,q[1])==jc-1:
                    below_black = (j==jc-1)
                    wgt*=np.exp(1j*theta) if below_black else np.exp(-1j*theta)
                K[bi[(i,j)],wi[q]]=wgt
    return K
def sector_distribution(L,H,M=64):
    th=2*np.pi*np.arange(M)/M; vals=[]
    s0,l0=np.linalg.slogdet(kast(L,H,0.0))
    for t in th:
        s,l=np.linalg.slogdet(kast(L,H,t)); vals.append(s/s0*np.exp(l-l0))
    vals=np.array(vals)
    # vals(theta) = sum_k p_k e^{i k theta}
    p=np.fft.fft(vals)/M      # p[k] for e^{-ik}... handle both
    ks=np.fft.fftfreq(M,1/M).astype(int)
    return {int(k):p[i] for i,k in enumerate(ks)}, l0
if __name__=='__main__':
    L,H=int(sys.argv[1]),int(sys.argv[2])
    P,l0=sector_distribution(L,H)
    print('L=%d hole=%d  ln(#tilings) = %.3f'%(L,H,l0))
    for k in sorted(P):
        if abs(P[k])>1e-6: print('  k=%+d  p=%.5f  (imag %.1e)'%(k,P[k].real,P[k].imag))
