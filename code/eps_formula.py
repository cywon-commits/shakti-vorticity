from fractions import Fraction as F
from sympy import factorint
vals={1:F(3,4),2:F(3*5**3*7,2**12),3:F(3**3*5**3*7**5*11,2**30),4:F(3**8*7**5*11**5*13**3,2**56),5:F(3**8*7*11**9*13**7*17**3*19,2**90),6:F(3**3*11**9*13**11*17**7*19**5*23,2**132),
7:F(2889099907996177643819591709751114247675619579860671875,6129982163463555433433388108601236734474956488734408704),
8:F(805490769012012267435003272840561406287076294522308802156245998212890625,1766847064778384329583297500742918515827483896875618958121606201292619776)}
if __name__=='__main__':
    for m,v in vals.items(): print(m,factorint(v.numerator),factorint(v.denominator))
    print('ratios')
    for m in range(1,8):
        r=vals[m+1]/vals[m]; print(m,'->',m+1,factorint(r.numerator),'/',factorint(r.denominator))
    print('second ratios')
    for m in range(1,7):
        r=(vals[m+2]/vals[m+1])/(vals[m+1]/vals[m]); print(m,factorint(r.numerator),'/',factorint(r.denominator))
