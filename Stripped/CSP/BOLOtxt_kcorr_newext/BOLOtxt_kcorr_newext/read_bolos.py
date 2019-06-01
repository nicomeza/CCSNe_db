import numpy as np 
import matplotlib.pyplot as pl
from astro_utils import * 
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
from scipy.optimize import newton,brentq


def Ni_K(Lpeak,tpeak,beta,Ni_min=0.01,Ni_max = 1.0):
    '''F =  lambda x : L_peak(x,Q_t,t_peak=tpeak,beta=beta)[0]-Lpeak                                                                                                   
    return brentq(F,Ni_min,Ni_max)'''
    F =  lambda x : L_peak(x,Q_t,t_peak=tpeak,beta=beta)[0]-Lpeak
    return brentq(F,Ni_min,Ni_max)

def L_peak(M_ni,L_heat,t_peak,beta=1.0):

    I = integrate.quad(lambda x: x*L_heat(x,M_ni), 0,beta*t_peak)
    return 2 * np.asarray(I) / (beta*t_peak)**2.



M_ni = lambda t_r,L_bol : L_bol/Q_t(t_r,1)

colormap = pl.cm.spectral

bolos = [l.split()[0] for l in open('Bolo.txt','r')]
filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(bolos))]
pl.gca().set_color_cycle(filters_colors)

for bolo in bolos:
    t,L = np.genfromtxt(bolo).T
    clean = np.where(L<1e43)[0]
    t,L = t[clean],L[clean]
    pl.plot(t,np.log10(L),marker='o',linestyle='None',label=bolo.split('_')[1])
    
pl.legend(loc='best',prop={'size':10},ncol=2)
pl.savefig("All_bolos.png")
pl.close()
#pl.show()

for bolo in bolos:
    sn_name = bolo.split('_')[1]
    print sn_name

    t,L = np.genfromtxt(bolo).T
    clean = np.where(L<1e43)[0]
    t,L = t[clean],L[clean]
    pl.plot(t,np.log10(L),marker='o',linestyle='None',label=sn_name,color='k')
    
    try:
     
        where_peak = np.where(np.logical_and(t>10.,t<45.))[0]
        xs,ys = t[where_peak],(L[where_peak])*1e-42
    #k0 = smooth.NonParamRegression(xs, ys, method=npr_methods.SpatialAverage())                                                                           
        k0 = smooth.NonParamRegression(xs, ys, method = npr_methods.LocalPolynomialKernel(q=4))
        k0_full = smooth.NonParamRegression(t,(L)*1e-42, method = npr_methods.LocalPolynomialKernel(q=7))
        k0.fit()
        k0_full.fit()
        grid = np.arange(np.min(xs),np.max(xs),0.1)
        grid_full = np.arange(np.min(t),np.max(t),0.1)
        y_grid = np.log10(k0(grid))+42
        y_grid = y_grid[~np.isnan(y_grid)]
        y_grid_full = np.log10(k0_full(grid_full))+42
        y_grid_full = y_grid_full[~np.isnan(y_grid_full)]
        where_peak = np.argmax(y_grid)
        tp,Lp = grid[where_peak],y_grid[where_peak]
        M_p = M_ni(tp,10**Lp)
        Ni_Khatami = Ni_K(10**Lp,tp,beta=4/3.)            
        ni_t = np.arange(5,np.max(t)+5,0.01)
            
        where_tail = np.where(t>tp)[0]
        x_tail,y_tail = t[where_tail][len(where_tail)-4:len(where_tail)],np.log10(L[where_tail][len(where_tail)-4:len(where_tail)])
        Ni_tail,Ni_tail_err = curve_fit(log_Q_t,x_tail,y_tail)
        print tp,Lp,x_tail,Ni_tail
        pl.plot(ni_t,np.log10(Q_t(ni_t,M_p)),linestyle='-.',color='k')
        pl.plot(ni_t,np.log10(Q_t(ni_t,Ni_tail)),linestyle='-.',color='k')
        pl.annotate(r"$\Delta ^{56}Ni = %2.2f$"%(M_p-Ni_tail),(x_tail[-1]-20,log_Q_t(x_tail[-1]-20,M_p)+0.1),size=20)

        pl.ylim(np.min(np.log10(L))-0.5,Lp+0.5)
        pl.xlabel(r'$\mathrm{Time \ since \ estimated \ explosion}$',size=15)
        pl.ylabel(r'$\mathrm{\log{L_{bol}}}$',size=20)
        pl.legend()
        pl.savefig('%s_bol.png'%sn_name)
        pl.close()
    except:
        
        pl.savefig('%s_bol.png'%sn_name)
        pl.close()
