from vort_rule import *
# coupling between the two dominoes of every 2x2 block (both fillings), anywhere in the lattice
C=Counter()
for x in range(-4,4):
    for y in range(-4,4):
        for fill in ('h','v'):
            if fill=='h': d1,d2=((x,y),(x+1,y)),((x,y+1),(x+1,y+1))
            else: d1,d2=((x,y),(x,y+1)),((x+1,y),(x+1,y+1))
            A,B=info2(d1),info2(d2); js,common,kinds=coupling(A,B)
            C[(fill,tuple(sorted(js)))]+=1
print(C)
