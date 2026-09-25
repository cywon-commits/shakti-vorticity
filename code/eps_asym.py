import sympy as sp, mpmath as mp
m=sp.symbols('m',positive=True); lnA=sp.Symbol('lnA')
def lnH(M): return (M**2/2+M/2+sp.Rational(1,12))*sp.log(M)-M**2/4+lnA        # hyperfactorial, + O(1/M^2)
def lnF(M): return M*sp.log(M)-M+sp.log(2*sp.pi*M)/2+1/(12*M)                  # Stirling
def S1(N): return lnH(2*N)-N*(N+1)*sp.log(2)-2*lnH(N)
def S0(N): return lnF(2*N)-N*sp.log(2)-lnF(N)
expr=4*S1(m)-S1(2*m)-8*m*S0(m)+4*m*S0(2*m)-2*m*(2*m-1)*sp.log(2)
expr=sp.expand(sp.expand_log(expr,force=True))
ser=sp.simplify(sp.series(expr.subs(m,1/sp.Symbol('t')),sp.Symbol('t'),0,1).removeO())
print('asymptotic ln|eps_m| =',sp.simplify(sp.collect(sp.expand(ser),sp.log(sp.Symbol('t')))))
# exact numeric check
mp.mp.dps=30
def exact_ln(M):
    s=mp.mpf(0)
    for j in range(1,4*M,2): s+=min(4*M-j,3*j-4*M)*mp.log(j)
    return s-2*M*(2*M-1)*mp.log(2)
Aval=mp.glaisher
for M in (8,100,1000,4000):
    e=exact_ln(M); print(M,'exact ln|eps| + (1/4)ln m =',mp.nstr(e+mp.log(M)/4,15))
