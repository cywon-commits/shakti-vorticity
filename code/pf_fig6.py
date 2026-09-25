import pickle, glob, math, numpy as np
from confrad import conformal_radius
from pubstyle import *
import matplotlib.ticker
def amp(D):
    W,H=D['W'],D['H']; V=D['vmean']; acc={}; cnt=0; r=min(W,H)/4
    for i in range(W):
        for j in range(H):
            if math.hypot(i+.5-W/2,j+.5-H/2)<=r:
                cnt+=1
                for q in [(1,1),(1,-1),(-1,1),(-1,-1)]:
                    acc[q]=acc.get(q,0)+V[i,j]*np.exp(-1j*math.pi/2*(q[0]*i+q[1]*j))
    return math.sqrt(sum(abs(a/cnt)**2 for a in acc.values()))
cr={}
rows=[]
for f in glob.glob('rect/r_*.pkl'):
    D=pickle.load(open(f,'rb')); W,H=D['W'],D['H']; a=W/H
    if a not in cr: cr[a]=conformal_radius(a,M=80)
    rows.append(dict(W=W,H=H,a=a,A=amp(D),rho=H*cr[a]))
rows.sort(key=lambda r:(r['a'],r['H']))
for r in rows: print('%3dx%-3d aspect %.2f  A=%.3f  H*rho(a)=%.2f'%(r['W'],r['H'],r['a'],r['A'],r['rho']))
x=np.log([r['rho'] for r in rows]); y=np.log([r['A'] for r in rows]); s,c=np.polyfit(x,y,1)
print('fit over all rectangles: A = %.3f * (H rho)^(%.3f)'%(math.exp(c),s))
xs=np.log([r['H'] for r in rows]); s2=np.polyfit(xs,y,1)[0]
print('naive fit A vs short side H only: exponent %.3f; residual std with rho: %.4f, with H: %.4f'%(s2,np.std(y-(s*x+c)),np.std(y-np.polyval(np.polyfit(xs,y,1),xs))))
fig,ax=plt.subplots(1,2,figsize=(COL,1.75),gridspec_kw=dict(wspace=0.45,left=0.14,right=0.99,top=0.88,bottom=0.23))
for k,a in enumerate((1,2,4)):
    R=[r for r in rows if r['a']==a]
    ax[0].plot([r['H'] for r in R],[r['A'] for r in R],MK[k],color=CAT[k],ms=3.5,label=f'$W/H={a:g}$')
HH=np.linspace(8,32,50)
for k,a in enumerate((1,2,4)): ax[0].plot(HH,math.exp(c)*(HH*cr[float(a)])**s,'-',color=CAT[k],lw=.6)
R=[r for r in rows if r['H']==16]; R.sort(key=lambda r:r['a'])
ax[1].plot([r['a'] for r in R],[r['A'] for r in R],'o',color=CAT[0],ms=3.5,label='MC, $H=16$')
aa=np.linspace(1,8,40); A1=[r['A'] for r in R if r['a']==1][0]
ax[1].plot(aa,[A1*(conformal_radius(a_,M=60)/cr[1.0])**-0.25 for a_ in aa],'k--',lw=.8,label='$\\rho^{-1/4}$ law')
ax[0].set_xscale('log'); ax[0].set_yscale('log'); ax[0].set_xlabel('short side $H$'); ax[0].set_ylabel('central amplitude')
ax[0].set_xticks([8,16,32]); ax[0].set_xticklabels(['8','16','32']); ax[0].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[0].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[0].set_yticks([.4,.5]); ax[0].set_yticklabels(['0.4','0.5'])
ax[0].legend(handletextpad=.2,borderaxespad=.2,loc='lower left'); ax[1].legend(loc='upper right',borderaxespad=.2)
ax[1].set_xlabel('aspect ratio $W/H$'); ax[1].set_ylabel('central amplitude')
ax[0].set_title('(a)',loc='left',pad=2); ax[1].set_title('(b)',loc='left',pad=2)
save(fig,'fig6_shape')
