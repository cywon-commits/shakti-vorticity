"""Phase of Psi vs flat (basketweave) domino arrangements, chessboard view (exploratory, not for the paper)."""
import sys, math, cmath, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, RegularPolygon, Rectangle
from pcolor import colour, wheel
from vform2 import S_TAB
from flat_phase import analyse, black
plt.rcParams.update({'font.family':'DejaVu Sans','mathtext.fontset':'dejavusans','font.size':8})
OR,BL,ORl,BLl='#e8743a','#2f6db5','#f6c6a8','#aac6ea'
CHOICE={0:(1,1,False,4),1:(0,1,False,0),2:(0,0,False,0),3:(1,0,False,0),4:(1,1,False,0),5:(0,1,False,4),6:(0,0,False,4),7:(1,0,False,4)}
fig=plt.figure(figsize=(7.0,7.6)); R0=0.36; SZ=0.205; Wn=4; YC=0.515; K=7.0/7.6
w=fig.add_axes([0.5-0.115,YC-0.115*K,0.23,0.23*K],projection='polar'); wheel(w,fs=7); w.set_xticklabels([])
for k in range(8):
    ph=k*math.pi/4
    x0,y0=0.5+0.165*math.cos(ph),YC+0.165*K*math.sin(ph); x1,y1=0.5+0.24*math.cos(ph),YC+0.24*K*math.sin(ph)
    fig.add_artist(plt.Line2D([x0,x1],[y0,y1],transform=fig.transFigure,color='k',lw=2.4 if k%2 else 1.0,solid_capstyle='butt'))
    al,be,par,off=CHOICE[k]; psi,hb,P,h=analyse(al,be,par,off,W=12)
    assert abs(((cmath.phase(psi)-ph+math.pi)%(2*math.pi))-math.pi)<1e-6
    cx,cy=0.5+R0*math.cos(ph),YC+R0*K*math.sin(ph)
    a=fig.add_axes([cx-SZ/2,cy-SZ*K/2,SZ,SZ*K]); a.set_xlim(-Wn-1.3,Wn+1.3); a.set_ylim(-Wn-1.3,Wn+1.3); a.set_aspect('equal'); a.axis('off')
    a.add_patch(RegularPolygon((0,0),8,radius=Wn+1.25,orientation=math.pi/8,fc=colour(0.30*cmath.exp(1j*ph)),ec='none',zorder=0))
    seen=set()
    for x in range(-Wn,Wn):
        for y in range(-Wn,Wn):
            d=tuple(sorted(((x,y),P((x,y)))))
            if d in seen or any(abs(c[0]+.5)>Wn or abs(c[1]+.5)>Wn for c in d): continue
            seen.add(d); p=d[0]; typ='H' if d[1][0]!=d[0][0] else 'V'
            v=S_TAB[(typ,p[0]%4,p[1]%4)]*(1 if (h[p]//4)%2==0 else -1)
            for c in d: a.add_patch(Rectangle(c,1,1,fc=(OR if black(c) else ORl) if v>0 else (BL if black(c) else BLl),ec='white',lw=.4,zorder=2))
            xs=[c[0] for c in d]; ys=[c[1] for c in d]
            a.add_patch(Rectangle((min(xs)+.06,min(ys)+.06),max(xs)-min(xs)+.88,max(ys)-min(ys)+.88,fill=False,ec='k',lw=.9,zorder=3))
    a.plot(0,0,'o',ms=5,mfc='white',mec='k',mew=1.0,zorder=5)
    lx,ly=0.5+0.143*math.cos(ph),YC+0.143*K*math.sin(ph)
    fig.text(lx,ly,f'{k}\n$\\frac{{{k}\\pi}}{{4}}$' if False else f'{k}',ha='center',va='center',fontsize=8.5,fontweight='bold')
    top=(k==2); a.text(0,(Wn+1.45) if top else (-Wn-1.45),f'$\\langle h\\rangle\\equiv{k}$, blocks $({al},{be})$',ha='center',va=('bottom' if top else 'top'),fontsize=7)
fig.text(0.5,0.995,'arg$\\,\\Psi=\\pi\\langle h\\rangle/4$: phase colour and the flat (basketweave) domino arrangements',ha='center',va='top',fontsize=9.5,fontweight='bold')
fig.text(0.015,0.015,'dominoes (chessboard view): orange $v=+1$, blue $v=-1$;\ndark/light = black/white cell;  white dot = reference vertex (0,0)\n'
         'blocks $(a,b)$: lower-left cells of the parallel-domino 2$\\times$2 blocks mod 2',ha='left',va='bottom',fontsize=6.8)
fig.text(0.985,0.015,'thick ticks: odd $\\langle h\\rangle$ (pinwheel)\nthin ticks: even $\\langle h\\rangle$ (facing)\n'
         '$k$ and $k+4$: same arrangement, all $v$ reversed',ha='right',va='bottom',fontsize=6.8)
fig.savefig('outputs/phase_wheel_chess.png',dpi=200); fig.savefig('outputs/phase_wheel_chess.pdf'); print('ok')
