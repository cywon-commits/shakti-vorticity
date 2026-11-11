"""Chessboard regions that correspond to axis-aligned rectangular Shakti samples."""
import io, contextlib, pickle, sys
with contextlib.redirect_stdout(io.StringIO()):
    from core import rect, pcen
from exact_v import Region, start_vertex, black
R0=70
CELLS={}
for i in range(-R0,R0):
    for j in range(-R0,R0):
        r=tuple(sorted(rect((i,j)))); CELLS[(i,j)]=r
def shakti_rect(x0,y0,X,Y):
    """all chess cells whose plaquette lies inside [x0,x0+X) x [y0,y0+Y) (Shakti unit squares)"""
    return {c for c,r in CELLS.items() if all(x0<=p[0]<x0+X and y0<=p[1]<y0+Y for p in r)}
if __name__=='__main__':
    X=Y=int(sys.argv[1]) if len(sys.argv)>1 else 24
    for ox in range(4):
        for oy in range(4):
            cells=shakti_rect(-X//2+ox,-Y//2+oy,X,Y)
            nb=sum(black(c) for c in cells); nw=len(cells)-nb
            try:
                Rg=Region(cells); ok='tileable, log Z=%.2f'%Rg.ld0
            except Exception as e: ok='NOT tileable (%s)'%type(e).__name__
            print(ox,oy,len(cells),'black',nb,'white',nw,ok,flush=True)
