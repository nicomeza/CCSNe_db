from astro_utils import * 
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
from scipy.misc import derivative
from scipy.optimize import newton,brentq
import itertools 
from matplotlib.lines import Line2D

prop_cycle = pl.rcParams['axes.color_cycle']
plot_pairs = list(itertools.product(prop_cycle,Line2D.filled_markers))

sn_labels = ['sn','type','host','hostredshift','hostlumdist','sn_ebv','sn_host','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_test_total.dat',dtype={'names':sn_labels,'formats':sn_formats})

compare_nis = False
get_peak = True
plot = True
show = False
Khatami = False
check_t_0 = False
check_type = False
Tail = True

M_ni = lambda t_r,L_bol : L_bol/Q_t(t_r,1)

current_dir = os.getcwd()

colormap = pl.cm.spectral

def L_peak(M_ni,L_heat,t_peak,beta=1.0):
    
    I = integrate.quad(lambda x: x*L_heat(x,M_ni), 0,beta*t_peak)
    return 2 * np.asarray(I) / (beta*t_peak)**2.

def Ni_K(Lpeak,tpeak,beta,Ni_min=0.01,Ni_max = 1.0):
    '''F =  lambda x : L_peak(x,Q_t,t_peak=tpeak,beta=beta)[0]-Lpeak
    return brentq(F,Ni_min,Ni_max)'''
    F =  lambda x : L_peak(x,Q_t,t_peak=tpeak,beta=beta)[0]-Lpeak
    return brentq(F,Ni_min,Ni_max)

def Ni_K2(Lpeak,tpeak,L_heat,beta=1.0):
    
    I = 2. * integrate.quad(lambda x: x*L_heat(x,1.0), 0,beta*tpeak)[0]
    return Lpeak * beta**2. * tpeak**2. / I 

def Ni_K3(Lpeak,tpeak,beta=1.0):
    t_ni = 8.8        # niquel decay [days]                                                                                       
    t_co = 111.3      # Cobalt decay [days]  
    x = beta*tpeak
    x_ni = x/t_ni
    x_co = x/t_co
    e_ni = 6.395*1e43
    e_co = 1.362*1e43
    corcho1 = e_ni*t_ni*(t_ni*(1-np.exp(-x_ni)) - x*np.exp(-x_ni))
    corcho2 = e_co*t_co*(t_co*(1-np.exp(-x_co)) - x*np.exp(-x_co))
    return x**2 * Lpeak / 2. / (corcho1+corcho2)

def weird_Ni(t,Mni,t_m=1.0):
    
    Weird = Q_t(t,Mni) - t*Q_t_dot(t,Mni)/((t/t_m)**2. - 1)
    return Weird

def weird_Ni_der(t,Mni,t_m=1.0):
    
    t_ni = 8.8        # niquel decay [days]                                                                                          
    t_co = 111.3      # Cobalt decay [days]                                                                                   
    Q = Mni*(6.45*np.exp(-t/t_ni)+1.45*np.exp(-t/t_co))*1e43  # [erg/s] . See Nadyozhin 1994   
    Q_dot = -Mni*(6.45*np.exp(-t/t_ni)/t_ni+1.45*np.exp(-t/t_co)/t_co)*1e43
    x = t/t_m
    Weird = 2*x**2. * Q_dot/(2.*x**2. - 1)
    return Weird

########################################################

#band_string = ("BVRIYJH_local","BVRIYJH_local")
band_string = ("BVRIYJH_shift","BVRIYJH_shift")
#band_string = ("BVRI","BVRcI")

inflection = False
pl.close("all")

if get_peak:
    
    if plot:
        fig1 = pl.figure(figsize=(14,6))
        ax1 = fig1.add_subplot(111)
        fig3 = pl.figure()
        ax3 = fig3.add_subplot(111)
        np.random.shuffle(plot_pairs)

    Lfile = open('Ni56_%s_peak.dat'%band_string[0],'w')
    Lfile.write("###############%s#######################\n"%band_string[0])
    Lfile.write("# SN \t SN type \t tpeak \t Lpeak [erg/s] \t log(Lpeak) \t M_ni \t M_ni_Khatami \t M_ni_tail t_1/2(<tp) \t t_1/2(>tp) \t FWHM \t t_inflection \t L_dot(inflection)\n")
    for SN,z_SN,E_B_V,t_0,d_L,sn_type,t_nd,t_d in SN_DATA[['sn','sn_z','sn_ebv','t_0','hostlumdist','type','t_non_det','t_discov']]:
        
        band_string = ("BVRIYJH_local","BVRIYJH_local")
        band_string = ("BVRIYJH_total","BVRIYJH_total")
#        band_string = ("BVRIJH","BVRIJH")
        #band_string = ("BVRI","BVRcI")
        print "####### %s ########## \n"%SN

        if plot:
            fig2 = pl.figure()
            ax2 = fig2.add_subplot(111)
            
        try: 
            os.chdir("%s/"%SN)
            t_half_1,t_half_2,der_at_inf,der_at_inf_2,Ni_tail,FWHM,Ni_tail_= None,None,None,None,None,None,None
            Tail = True
            if t_0<1000.:
                t_0 = 0.5*(t_d+t_nd)
            try:
                
                try:
                    #print band_string
                    
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[0])).T
                    where_nan = np.intersect1d(np.where(~np.isnan(L))[0],np.unique(t,return_index=True)[1])
                    t,L,logL,Mni = t[where_nan],L[where_nan],logL[where_nan],Mni[where_nan]
                    band_string = band_string[0]
                    print '%s_Lbol_%s.dat'%(SN,band_string)
                except:
                    #print band_string
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[1])).T
                    where_nan = np.intersect1d(np.where(~np.isnan(L))[0],np.unique(t,return_index=True)[1])
                    t,L,logL,Mni = t[where_nan],L[where_nan],logL[where_nan],Mni[where_nan] 
                    band_string = band_string[1]
                    print '%s_Lbol_%s.dat'%(SN,band_string)

                t = t/(1+z_SN)
                xs,ys = t[np.argsort(t)],logL[np.argsort(t)]
                t,logL = xs,ys
                dt_0 = (t_0-(t_d-1))
                if dt_0 > 0 :
                    dt_0 = 0.

                print sn_type
                if sn_type=="IIb":
                    where_peak = np.where(np.logical_and(xs>10.,xs<55.))[0]
                else:
                    where_peak = np.where(xs<55.)[0]

                if sn_type in ["Ib","Ic","Ibc","Ic_GRB","Ic_BL"]:
                    beta_SN = 9/8.
                else:
                    beta_SN = 0.82

                xs,ys = xs[where_peak],(10**ys[where_peak])*1e-42
                #k0 = smooth.NonParamRegression(xs, ys, method=npr_methods.SpatialAverage())
                k0 = smooth.NonParamRegression(xs, ys, method = npr_methods.LocalPolynomialKernel(q=2))
                if len(t)>9:
                    k0_full = smooth.NonParamRegression(t,(10**logL)*1e-42, method = npr_methods.LocalPolynomialKernel(q=7))
                else:
                    k0_full = smooth.NonParamRegression(t,(10**logL)*1e-42, method = npr_methods.LocalPolynomialKernel(q=3))
                k0.fit()
                k0_full.fit()
                k0_res = k0(xs)-ys
                k0_res_full = k0_full(t)-(10**logL)*1e-42
                #pl.plot(xs,k0_res)
                #pl.show()
                
                #print k0_res,k0_res_full
                
                grid = np.arange(np.min(xs),np.max(xs),0.2)
                grid_full = np.arange(np.min(t),np.max(t),0.2)
                y_grid = np.log10(k0(grid))+42
                y_grid,grid = y_grid[~np.isnan(y_grid)],grid[~np.isnan(y_grid)]
                y_grid_full = np.log10(k0_full(grid_full))+42
                y_grid_full,grid_full = y_grid_full[~np.isnan(y_grid_full)],grid_full[~np.isnan(y_grid_full)]
                where_peak = np.argmax(y_grid)
                tp,Lp = grid[where_peak],y_grid[where_peak]
                M_p = M_ni(tp,10**Lp)
                Ni_Khatami = Ni_K(10**Lp,tp,beta=beta_SN)
                M_p_ = M_ni(tp+dt_0,10**Lp)
                Ni_Khatami_ = Ni_K(10**Lp,tp+dt_0,beta=beta_SN)
                ni_t = np.arange(np.argmin([tp-20,0]),np.max(t)+5,0.01)
                print tp,Lp,M_p,Ni_Khatami,M_p_,Ni_Khatami_
                after_p,before_p = np.where(grid>tp+2)[0],np.where(grid<=tp-2)[0]
                after_p_full,before_p_full = np.where(grid_full>tp+2)[0],np.where(grid_full<=tp-2)[0]
                order = 5

                if len(before_p):

                    print "Before peak :\n"

                    try:
                        
                        f = lambda x : k0_full(x) - 10**(Lp-42)/2.
                        zero = brentq(f,np.min(t),tp-2)     
                        t_half_1 = tp - zero 
                        print "t_1/2 = %s"%t_half_1
                    except:
                        print "No t_1/2"

                    try:
                        der2 = derivative(k0_full, grid_full[before_p_full], dx=0.5, n=2, order=order)
                        inf_t = grid_full[before_p_full][np.argmin(np.abs(der2))]
                        der_at_inf =  1e42*derivative(k0, inf_t, dx=0.5, n=2, order=order)
                        print (der_at_inf - Q_t_dot(inf_t,M_p))/Q_t_dot(inf_t,M_p)
                    
                        if plot and inflection:
                            for t_m in [5,10,15,20,25,30]:
                                #print t_m,(weird_Ni_der(inf_t,M_p,t_m=t_m) - der_at_inf)/der_at_inf
                                pl.axvline(inf_t,color='k',linestyle='--')
                                pl.annotate("Inflection",(inf_t-2.,np.min(logL)+0.1),size=8,rotation='vertical')
                    except:
                        print "No inflection"

                if len(after_p):
                    print "After peak :\n"

                    try:
                        
                        f = lambda x : k0_full(x) - 10**(Lp-42)/2.
                        zero = brentq(f,tp+2,np.max(t))     
                        t_half_2 = zero - tp
                        print "t_1/2 = %s"%t_half_2
                    except:
                        print "No t_1/2"

                    try:
                        der2 = derivative(k0, grid[after_p], dx=0.5, n=2, order=order)
                        inf_t2 = grid[after_p][np.argmin(np.abs(der2))]
                        der_at_inf_2 =  1e42 * derivative(k0, inf_t2, dx=0.5, n=2, order=order)
                        print (der_at_inf_2 - Q_t_dot(inf_t2,M_p))/Q_t_dot(inf_t2,M_p)
                    
                        for t_m in [5,10,15,20,25,30]:
                            print t_m,(weird_Ni_der(inf_t2,M_p,t_m=t_m) - der_at_inf_2)/der_at_inf_2                    
                            if plot and inflection:
                                pl.axvline(inf_t2,color='k',linestyle='--')
                                pl.annotate("Inflection",(inf_t2-2.,np.min(logL)+0.1),size=8,rotation='vertical')
                    except:
                        print "No inflection"
                
                try:
                    if plot:
                        if inflection:
                            for t_m,color in zip([5,10,20,30],['g','r','b','orange']):
                                ax2.plot(ni_t,np.log10(weird_Ni(ni_t,M_p,t_m=t_m)),linestyle='--',color=color,label='t_m = %s'%t_m)                    
                except:
                    print "No inflection"
                    
                try:
                    FWHM = t_half_2+t_half_1
                    print "FWHM: %s"%FWHM
                except:
                    FWHM = 99


                if Tail:
                    try: 
                        Tail = False
                        if t_half_2>0.:
                            where_tail = np.where(t>np.min([tp+2*t_half_2,tp+25.]))[0]
                        else:
                            where_tail = np.where(t>tp+25)[0]
                        if len(where_tail)>2:
                            x_tail,y_tail = t[where_tail][len(where_tail)-3:len(where_tail)],logL[where_tail][len(where_tail)-3:len(where_tail)]
                            Ni_tail,Ni_tail_err = curve_fit(log_Q_t,x_tail,y_tail)
                            Ni_tail,Ni_tail_err = Ni_tail[0],Ni_tail_err[0][0]
                            Tail = True
                            print "Tail fit : %s +- %s"%(Ni_tail,Ni_tail_err)
                            Ni_tail_,Ni_tail_err_ = curve_fit(log_Q_t,x_tail+dt_0,y_tail)
                            Ni_tail_,Ni_tail_err_ = Ni_tail_[0],Ni_tail_err_[0][0]
                            
                        else:
                            print "Tail fit failed"
                            Tail = False
                    except:
                        print "Tail fit failed"
                

                if plot:

                    color,marker = plot_pairs[np.where(SN_DATA['sn']==SN)[0][0]]
                    ax2.plot(grid_full,y_grid_full,linestyle='--',color='g',linewidth=4)
                    #ax2.plot(grid,y_grid,color='b')
                    ax2.plot(t,logL,marker='o',color='k',linestyle='None',label=r"$\mathrm{%s}$"%SN,alpha=0.4,markersize=10)
                    if sn_type=="Ic_BL" or sn_type=="Ic_GRB":
                        ax1.plot(t,logL,marker=marker,color=color,linestyle='None',\
                                     label=r"$\mathrm{%s (%s-%s)}$"%(SN[0:2]+SN[4:],sn_type.split('_')[0],sn_type.split('_')[1]),alpha=0.6)
                    else:
                        ax1.plot(t,logL,marker=marker,color=color,linestyle='None',label=r"$\mathrm{%s (%s)}$"%('SN'+SN[4:],sn_type),alpha=0.4,markersize=9)

                    ax2.plot(ni_t,np.log10(Q_t(ni_t,M_p)),linestyle='-.',color='k',linewidth=3) # Arnett's rule
                    ax2.plot(ni_t,np.log10(Q_t(ni_t,Ni_Khatami)),linestyle='--',color='k',linewidth=3) # Khatami rule
                    
                    ax1.plot(ni_t,np.log10(Q_t(ni_t,0.26)),linestyle='--',color='k',linewidth=2,alpha=0.7) 
                    
                    if Tail:

                        ax2.plot(ni_t,np.log10(Q_t(ni_t,Ni_tail)),linestyle=':',color='k',linewidth=3)  # Tail value
                        ax2.annotate(r"$^{56}\mathrm{Ni} \ L_{peak} = \ %2.2f$"%(M_p),(x_tail[-1]*0.7,log_Q_t(np.max([x_tail[-1]*0.7,20]),M_p)+0.2),size=22)
                        ax2.annotate(r"$^{56}\mathrm{Ni} \ L_{tail} \ = \ %2.2f$"%(Ni_tail),(x_tail[-1]*0.7,log_Q_t(np.max([x_tail[-1]*0.7,20]),M_p)+0.12),size=22)
                        ax2.annotate(r"$^{56}\mathrm{Ni} \ L_{k&k} \ = \ %2.2f$"%(Ni_Khatami),(x_tail[-1]*0.7,log_Q_t(np.max([x_tail[-1]*0.7,20]),M_p)+0.04),size=22)
                        #ax3.plot(t,L/Q_t(t,Ni_tail),linestyle='--',color='k',label=SN)
                        #L_dot = 1e42*derivative(k0_full,grid_full, dx=0.5, n=1, order=order)
                        k0_tail = smooth.NonParamRegression(t,logL, method = npr_methods.LocalPolynomialKernel(q=2))
                        k0_tail.fit()
                        L_dot = np.gradient(10**k0_tail(t),t)
                        Q_dot = Q_t_dot(t,Ni_tail)
                        ax3.plot(t-tp,(L_dot-Q_dot)/(-Q_dot),linestyle='--',color=color,label=r"$\mathrm{%s}$"%SN,marker=marker,alpha=0.5,markersize=9)
                    else:
                        ax2.annotate(r"$^{56}\mathrm{Ni} \ L_{peak} = \ %2.2f$"%(M_p),(t[-1]*0.7,log_Q_t(np.max([t[-1]*0.7,20]),M_p)+0.2),size=22)
                    #pl.plot(ni_t,np.log10(Q_t(ni_t,Ni_Khatami)),linestyle='-.',color='k')
                    #pl.axhline(Lp-np.log10(2),linestyle='--',color='k')
                    #pl.annotate("Half max",(np.max(t)-5,Lp-np.log10(2)+0.05),size=8)
                    

                if plot:

                    ax2.set_xlabel(r'$\mathrm{Time \ since \ estimated \ explosion}$',size=22)
                    ax2.set_ylabel(r'$\mathrm{\log{L_{bol}(%s)}}$'%(band_string.split('_')[0]),size=22)
                    ax2.legend(prop={'size':17})
                    ax2.set_ylim(np.min(logL)-0.1,Lp+0.2)
                    ax2.tick_params(labelsize=17)
                    fig2.subplots_adjust(top=0.98,bottom=0.11,right=0.97)
                    fig2.savefig("../%s_%s_Lbol.png"%(SN,band_string),dpi=400)
                    
                    if show:
                        pl.show()
                        
                        
                    fig4 = pl.figure()
                    ax4 = fig4.add_subplot(111)
                    ax4.plot(t,Mni,marker='o',color='y',linestyle='None')
                    ax4.plot(grid,M_ni(grid,10**y_grid),color='y',linestyle='--')
                    ax4.set_xlabel('t-t_0')
                    ax4.set_ylabel('M(Ni56)')
                    ax4.axvline(tp,color='k')
                    ax4.axvline(inf_t2,color='k')
                    ax4.axhline(M_p,color='k')
                    ax4.set_xlim(np.max([0,np.min(t)-5]),np.max(t)+5)
                    ax4.tick_params(labelsize=17)
                    fig4.savefig('%s_%s_Ni.png'%(SN,band_string))
                    if show:
                        pl.show()
                    else:
                        pl.close("all")
                        
                s = "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n"\
                    %(SN,sn_type,tp,10**Lp,Lp,M_p,Ni_Khatami,Ni_tail,t_half_1,t_half_2,FWHM,inf_t2,der_at_inf_2*1e-38)
                print s 
                s = "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n"\
                    %(SN,sn_type,tp+dt_0,10**Lp,Lp,M_p_,Ni_Khatami_,Ni_tail_,t_half_1,t_half_2,FWHM,inf_t2,der_at_inf_2*1e-38)
                print s 
                Lfile.write(s)
                os.chdir(current_dir)
            except:
                print "No Lbol file"
                os.chdir(current_dir)
        except:
            print "No directory %s/"%SN
    
    if plot:

        ax1.set_xlabel(r'$\mathrm{Time \ since \ estimated \ explosion}$',size=22)
        ax1.set_ylabel(r'$\mathrm{\log{L_{bol}(%s)}}$'%(band_string.split('_')[0]),size=22)
        ax1.set_ylim(41,42.75)
        ax1.set_xlim(0,90)
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.6, box.height])
        ax1.legend(ncol=2,bbox_to_anchor=(1.61,1.),loc='upper right',prop={'size':8.45})
        pl.setp(fig1.gca().get_legend().get_texts(), fontsize='16') 
        ax1.tick_params(labelsize=17)
        fig1.subplots_adjust(left=0.08)
        fig1.subplots_adjust(right=0.65)
        fig1.subplots_adjust(top=0.98,bottom=0.11)
        
        fig1.savefig("All_bolos.pdf",format='pdf',dpi=1000)

        ax3.axhline(0.0,linestyle='--',color='k')
        ax3.axvline(0.0,linestyle='--',color='k')
        ax3.set_xlabel(r'$\mathrm{Time \ since \ estimated \ peak}$',size=22)
        ax3.set_ylabel(r'$\mathrm{(\dot{L}_{bol}(%s)-\dot{Q_t})}/|\dot{Q_t}|}$'%(band_string.split('_')[0]),size=22)
        box = ax3.get_position()
        ax3.set_position([box.x0, box.y0, box.width * 0.96, box.height])
        ax3.legend(ncol=1,bbox_to_anchor=(1.4,1.065),loc='upper right',prop={'size':11})
        ax3.tick_params(labelsize=17)
        fig3.subplots_adjust(left=0.12)
        fig3.subplots_adjust(right=0.75)
        fig3.subplots_adjust(top=0.95)
        pl.setp(fig3.gca().get_legend().get_texts(), fontsize='16') 
        fig3.savefig("All_decay.pdf",format='pdf',dpi=1000)

    Lfile.close()

