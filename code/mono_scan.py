import pickle, random
from mono_tile import *
L=22; board={(i,j) for i in range(L) for j in range(L)}; rng=random.Random(5)
a=(6,11); res={'mono':[], 'O':[]}
for dx in range(1,12):
    b=(a[0]+dx,a[1]+(0 if dx%2 else 1))
    d=abs(b[0]-a[0])+abs(b[1]-a[1])
    fs=[]
    for t in range(8):
        tiles=random_tiling(board-{a,b},rng); J,_=couplings(tiles); fs.append(min_frustration(len(tiles),J)[0])
    res['mono'].append((d,fs)); print('mono',d,fs,flush=True)
# O pairs: O at a (cells a, partner(a)) and O at b
for dx in range(2,12,2):
    b=(a[0]+dx,a[1]+1)
    Oa=tuple(sorted((a,partner(a)))); Ob=tuple(sorted((b,partner(b))))
    if set(Oa)&set(Ob): continue
    fs=[]
    for t in range(8):
        try: tiles=random_tiling(board-set(Oa)-set(Ob),rng)
        except AssertionError: break
        tiles+= [('O',)+Oa,('O',)+Ob]; J,conf=couplings(tiles); fs.append(min_frustration(len(tiles),J)[0])
    res['O'].append((dx,fs)); print('O',dx,fs,flush=True)
pickle.dump(res,open('mono/tile_scan.pkl','wb'))
