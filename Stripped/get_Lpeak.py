from astro_utils import * 
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
from scipy.misc import derivative
from scipy.optimize import newton,brentq
import itertools 
from matplotlib.lines import Line2D

prop_cycle = pl.rcParams['axes.color_cycle']
plot_pairs = list(itertools.product(prop_cycle,Line2D.filled_markers))

sn_labels = ['sn','type','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_test.dat',dtype={'names':sn_labels,'formats':sn_formats})

compare_nis = False
get_peak = True
plot = True
show = False
Khatami = False
check_t_0 = False
check_type = False
Tail = False

M_ni = lambda t_r,L_bol : L_bol/Q_t(t_r,1)

current_dir = os.getcwd()

colormap = pl.cm.spectral

def L_peak_explicit(x,Mni=0.1,t_peak=20):

    A,B,C,D,E,F = 1.0402,-2.4935,-84.23,-0.28336,-0.75684,86.72
    t_ni = 8.8        # niquel decay [days]                                                                                       
    t_co = 111.3      # Cobalt decay [days]  
    x_ni = x/t_ni
    x_co = x/t_co

    corcho = (A*x + B*np.exp(-x_ni) + C*np.exp(-x_co) + D*x*np.exp(-x_ni) + E*x*np.exp(-x_co) + F) * 1e12 /x**2
    return 2*Mni*corcho

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
band_string = ("BVRIYJH_total","BVRIYJH_total")
#band_string = ("BVRI","BVRcI")

inflection = False
pl.close("all")

if get_peak:
    
    if plot:
        fig1 = pl.figure()
        ax1 = fig1.add_subplot(111)
        np.random.shuffle(plot_pairs)

    Lfile = open('Ni56_%s_peak.dat'%band_string[0],'w')
    Lfile.write("###############%s#######################\n"%band_string[0])
    Lfile.write("# SN \t SN type \t tpeak \t Lpeak [erg/s] \t log(Lpeak) \t M_ni \t M_ni_Khatami \t M_ni_tail t_1/2(<tp) \t t_1/2(>tp) \t FWHM \t t_inflection \t L_dot(inflection)\n")
    for SN,z_SN,E_B_V,t_0,d_L,sn_type in SN_DATA[['sn','sn_z','sn_ebv','t_0','hostlumdist','type']]:
        
        band_string = ("BVRIYJH_local","BVRIYJH_local")
        band_string = ("BVRIYJH_total","BVRIYJH_total")
        #band_string = ("BVRI","BVRcI")
        print "####### %s ########## \n"%SN

        if plot:
            fig2 = pl.figure()
            ax2 = fig2.add_subplot(111)
    
        try: 
            os.chdir("%s/"%SN)
            t_half_1,t_half_2,der_at_inf,der_at_inf_2,Ni_tail= None,None,None,None,None
            Tail = True
            try:
                
                try:
                    print band_string
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[0])).T
                    where_nan = np.where(~np.isnan(L))[0]
                    t,L,logL,Mni = t[where_nan],L[where_nan],logL[where_nan],Mni[where_nan] 
                    band_string = band_string[0]
                except:
                    print band_string
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[1])).T
                    band_string = band_string[1]

                t = t/(1+z_SN)
                xs,ys = t[np.argsort(t)],logL[np.argsort(t)]
                t,logL = xs,ys
                print sn_type
                if sn_type=="IIb":
                    where_peak = np.where(np.logical_and(xs>10.,xs<45.))[0]
                else:
                    where_peak = np.where(xs<45.)[0]

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
                ni_t = np.arange(np.argmin([tp-20,0]),np.max(t)+5,0.01)
                print tp,Lp,M_p,Ni_Khatami

                if Tail:
                    try: 
                        Tail = False
                        where_tail = np.where(t>tp+25)[0]
                        if len(where_tail)>2:
                            x_tail,y_tail = t[where_tail][len(where_tail)-3:len(where_tail)],logL[where_tail][len(where_tail)-3:len(where_tail)]
                            Ni_tail,Ni_tail_err = curve_fit(log_Q_t,x_tail,y_tail)
                            Ni_tail,Ni_tail_err = Ni_tail[0],Ni_tail_err[0][0]
                            Tail = True
                            print "Tail fit : %s +- %s"%(Ni_tail,Ni_tail_err)
                        else:
                            print "Tail fit failed"
                            Tail = False
                    except:
                        print "Tail fit failed"
                

                if plot:

                    color,marker = plot_pairs[np.where(SN_DATA['sn']==SN)[0][0]]
                    ax2.plot(grid_full,y_grid_full,linestyle='--',color='g')
                    ax2.plot(grid,y_grid,color='b')
                    ax2.plot(t,logL,marker='o',color='k',linestyle='None',label=r"$\mathrm{%s}$"%SN,alpha=0.6)
                    ax1.plot(t,logL,marker=marker,color=color,linestyle='None',label=r"$\mathrm{%s}$"%SN,alpha=0.6)
                    ax2.plot(ni_t,np.log10(Q_t(ni_t,M_p)),linestyle='-.',color='k')
                    
                    if Tail:

                        ax2.plot(ni_t,np.log10(Q_t(ni_t,Ni_tail)),linestyle='-.',color='k')
                        ax2.annotate(r"$^{56}\mathrm{Ni} \ L_{peak} = \ %2.2f$"%(M_p),(x_tail[-1]*0.7,log_Q_t(np.max([x_tail[-1]*0.7,20]),M_p)+0.2),size=20)
                        ax2.annotate(r"$^{56}\mathrm{Ni} \ L_{tail} \ = \ %2.2f$"%(Ni_tail),(x_tail[-1]*0.7,log_Q_t(np.max([x_tail[-1]*0.7,20]),M_p)+0.12),size=20)
                    else:
                        ax2.annotate(r"$^{56}\mathrm{Ni} \ L_{peak} = \ %2.2f$"%(M_p),(t[-1]*0.7,log_Q_t(np.max([t[-1]*0.7,20]),M_p)+0.2),size=20)
                    #pl.plot(ni_t,np.log10(Q_t(ni_t,Ni_Khatami)),linestyle='-.',color='k')
                    #pl.axhline(Lp-np.log10(2),linestyle='--',color='k')
                    #pl.annotate("Half max",(np.max(t)-5,Lp-np.log10(2)+0.05),size=8)
                    
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

                if plot:

                    ax2.set_xlabel(r'$\mathrm{Time \ since \ estimated \ explosion}$',size=15)
                    ax2.set_ylabel(r'$\mathrm{\log{L_{bol}(%s)}}$'%(band_string),size=20)
                    ax2.legend(prop={'size':15})
                    ax2.set_ylim(np.min(logL)-0.1,Lp+0.2)
                    fig2.savefig("../%s_%s_Lbol.png"%(SN,band_string))
                    
                    if show:
                        pl.show()

                    pl.figure()
                    pl.plot(t,Mni,marker='o',color='y',linestyle='None')
                    pl.plot(grid,M_ni(grid,10**y_grid),color='y',linestyle='--')
                    pl.xlabel('t-t_0')
                    pl.ylabel('M(Ni56)')
                    pl.axvline(tp,color='k')
                    pl.axvline(inf_t2,color='k')
                    pl.axhline(M_p,color='k')
                    pl.xlim(np.max([0,np.min(t)-5]),np.max(t)+5)
                    pl.savefig('%s_%s_Ni.png'%(SN,band_string))
                    if show:
                        pl.show()
                    else:
                        pl.close("all")
                    
                s = "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n"\
                    %(SN,sn_type,tp,10**Lp,Lp,M_p,Ni_Khatami,Ni_tail,t_half_1,t_half_2,FWHM,inf_t2,der_at_inf_2*1e-38)
                print s 
                Lfile.write(s)
                os.chdir(current_dir)
            except:
                print "No Lbol file"
                os.chdir(current_dir)
        except:
            print "No directory %s/"%SN
    
    if plot:
        ax1.set_xlabel(r'$\mathrm{Time \ since \ estimated \ explosion}$',size=15)
        ax1.set_ylabel(r'$\mathrm{\log{L_{bol}(%s)}}$'%(band_string),size=20)
        ax1.legend(prop={'size':9},ncol=2,bbox_to_anchor=(1.15, 1))
        fig1.savefig("All_bolos.png")
    Lfile.close()


prentice_file = "Full_stats.dat"
my_file = "Ni56_BVRIYJH_peak.dat"

if compare_nis:

    pl.close("all")
    if prentice_file == "Full_stats.dat":
        nis = np.genfromtxt('Prentice/%s'%prentice_file,\
                                names=('SN','type','logLp','logLp_err1','logLp_err2','tp','tp_err','Mni','Mni_err1','Mni_err2'),\
                                dtype=('S15','S10','f8','f8','f8','f8','f8','f8','f8','f8'))
        compare_items = [('Mni'),('logLp'),('tp')]
    else:
        nis = np.genfromtxt('Prentice/%s'%prentice_file,\
                                names=('SN','type','logLp','logLp_err1','logLp_err2','Mni','Mni_err1','Mni_err2','tp','tp_err','trise','trise_err','tdecay','tdecay_err','FWHM','FWHM_err'),\
                                dtype=('S15','S10','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
        compare_items = [('Mni'),('logLp'),('tp'),('trise'),('tdecay'),('FWHM')]
    my_nis = np.genfromtxt(my_file,\
                               names=('SN','type','tp','Lp','logLp','Mni','Mni_K','Mni_tail','trise','tdecay','FWHM','t_inf','L_dot_inf'),\
                               dtype=('S15','S10','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

    i1,i2 = overlap(my_nis['SN'],nis['SN'])
    
    
    for item in compare_items:
        print item 
        pl.figure()
        no_99 = np.logical_and(nis[i2][item] != 99.,my_nis[i1][item] != 99.)
        no_nan = np.logical_and(~np.isnan(my_nis[i1][item]),no_99)
        if len(np.where(no_nan)[0])>0:

            if item == 'Mni' or item == "logLp":
                pl.errorbar(my_nis[i1][item][no_nan],nis[i2][item][no_nan],yerr = (nis[i2]['%s_err1'%item][no_nan],nis[i2]['%s_err2'%item][no_nan]),fmt='o')
            else:
                pl.errorbar(my_nis[i1][item][no_nan],nis[i2][item][no_nan],nis[i2]['%s_err'%item][no_nan],fmt='o')
            
            xmin,xmax = np.min([my_nis[i1][item][no_nan],nis[i2][item][no_nan]]),np.max([my_nis[i1][item][no_nan],nis[i2][item][no_nan]])
            x = np.arange(xmin-0.3*(xmax-xmin),0.3*(xmax-xmin)+xmax,0.01)
            print xmin,xmax
            pl.xlabel('My %s'%item)
            pl.ylabel('Prentice %s'%item)
            pl.plot(x,x)
            pl.xlim(xmin-0.3*(xmax-xmin),xmax+0.3*(xmax-xmin))
            pl.ylim(xmin-0.3*(xmax-xmin),xmax+0.3*(xmax-xmin))
            
            for sn,mni,mni2 in zip(my_nis[i1]['SN'][no_nan],my_nis[i1][item][no_nan],nis[i2][item][no_nan]):
                pl.annotate(sn,(mni+0.005,mni2+0.005),size=7)
            pl.savefig("%s_%s_%s_comp.png"%(item,my_file.split('_')[1],prentice_file.split('_')[0]))
        else:
            pl.close()
    if show:
        pl.show()
    else:
        pl.close()

    scatter = True
    hist = False
    for item1,item2 in [('Mni','Mni_K')]:
        if scatter:
            pl.figure()
            no_99 = np.logical_and(my_nis[item1] != 99.,my_nis[item2] != 99.)
            no_nan = np.logical_and(~np.isnan(my_nis[item1]),~np.isnan(my_nis[item2]))
            valid = np.where(np.logical_and(no_99,no_nan))[0]
            xmin,xmax = np.min(my_nis[item1][valid]),np.max(my_nis[item1][valid])
            ymin,ymax = np.min(my_nis[item2][valid]),np.max(my_nis[item2][valid])
            xymin,xymax = np.min([xmin,ymin]),np.max([xmax,ymax])
            x = np.arange(xymin-0.3*(xymax-xymin),0.3*(xymax-xymin)+xymax,0.01)
            pl.plot(x,x)
            pl.plot(my_nis[item1][valid],my_nis[item2][valid],'o')
            pl.xlim(xymin-0.3*(xymax-xymin),xymax+0.3*(xymax-xymin))
            pl.xlabel(item1)
            pl.ylabel(item2)

        if hist:
            
            pl.figure()
            no_99 = np.logical_and(my_nis[item1] != 99.,my_nis[item2] != 99.)
            no_nan = np.logical_and(~np.isnan(my_nis[item1]),~np.isnan(my_nis[item2]))
            valid = np.where(np.logical_and(no_99,no_nan))[0]
            pl.hist(my_nis[item1][valid],label=item1)
            pl.hist(my_nis[item2][valid],label=item2)
            pl.legend()

    pl.show()





if Khatami:    

    betas = np.arange(0.6,2.1,0.05)

    t_peak = 20
    Lps = np.arange(41.5,43,0.5)

    filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(Lps))]       
    pl.gca().set_color_cycle(filters_colors)
    
    
    for Lp in Lps:
        Ni_K_0 = Ni_K(10**Lp,t_peak,beta=1.0,Ni_min=0.005,Ni_max = 1.0)
        Ks = []
        
        for beta in betas:
            
            try:
                Ks.append(Ni_K(10**Lp,t_peak,beta=beta,Ni_min=0.005,Ni_max = 1.0))
            except:
                print "failed at %s"%beta
                Ks.append(0.0)
        pl.plot(betas,Ks,label=r'$log(L_{peak} = %2.2f)$'%Lp)


        pl.plot(0,0,label=r'$t_{peak} = %s$'%t_peak,linestyle='--',color='k')

    t_peak = 40
    filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(Lps))]       
    pl.gca().set_color_cycle(filters_colors)
    
    for Lp in Lps:

        Ni_K_0 = Ni_K(10**Lp,t_peak,beta=1.0,Ni_min=0.005,Ni_max = 1.0)
        Ks = []
        
        for beta in betas:
            
            try:
                Ks.append(Ni_K(10**Lp,t_peak,beta=beta,Ni_min=0.005,Ni_max = 1.0))
            except:
                print "failed at %s"%beta
                Ks.append(0.0)
        pl.plot(betas,Ks,linestyle='--')


    pl.plot(0,0,label=r'$t_{peak} = %s$'%t_peak,linestyle='--',color='k')
    pl.xlabel(r'$\beta$')
    pl.ylabel(r'$M(^{56}\mathrm{Ni})$')
    pl.legend(loc='best',prop={'size':10},ncol=2)
    pl.show()


    Lp = 42.5
    t_ps = np.arange(5,40,5)
    filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(t_ps))]       
    pl.gca().set_color_cycle(filters_colors)
    
    for tp in t_ps:
        
        Ks = []
        
        for beta in betas:
            
            try:
                Ks.append(Ni_K(10**Lp,tp,beta=beta,Ni_min=0.005,Ni_max = 1.0)/Ni_K(10**Lp,tp,beta=1.0,Ni_min=0.005,Ni_max = 1.0))
            except:
                print "failed at %s"%beta
                Ks.append(0.0)

        pl.plot(betas,Ks,linestyle='-',label=r'$t_{peak} = %2.1f$'%tp)
        
    #pl.plot(0,0,label=r'$log(L_{peak}) = %s$'%Lp,linestyle='-',color='k')
        
    for beta_SN,beta_label in [(9/8.,"Ibc"),(0.82,"IIb/pec"),(2.0,'Generic/Mixed'),(0.7,'Central/Helium')]:
        pl.axvline(beta_SN,linestyle='--',linewidth=4,color='k',alpha=0.7)
        pl.annotate(r"$\mathrm{%s}$"%(beta_label),(beta_SN-0.07,1.3),size=20,rotation=90)

    pl.xlabel(r'$\beta$',size=18)
    pl.ylabel(r'$M(^{56}\mathrm{Ni})/M(^{56}\mathrm{Ni})(\beta=1.0)$',size=18)
    pl.legend(loc='best',prop={'size':15},ncol=2)
    pl.show()

    
    for SN,z_SN,E_B_V,t_0,d_L,sn_type in SN_DATA[['sn','sn_z','sn_ebv','t_0','hostlumdist','type']]:
    
        band_string = ("BVRIYJH","BVRIYJH")
        #band_string = ("BVRI","BVRcI")
        print "####### %s ########## \n"%SN
        
        try: 
            os.chdir("%s/"%SN)
            t_half_1,t_half_2,der_at_inf,der_at_inf_2,Ni_tail= None,None,None,None,None
            Tail = True
            try:
                
                try:
                    print band_string
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[0])).T
                    where_nan = np.where(~np.isnan(L))[0]
                    t,L,logL,Mni = t[where_nan],L[where_nan],logL[where_nan],Mni[where_nan] 
                    band_string = band_string[0]
                except:
                    print band_string
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[1])).T
                    band_string = band_string[1]

                t = t/(1+z_SN)
                xs,ys = t[np.argsort(t)],logL[np.argsort(t)]
                t,logL = xs,ys
                print sn_type
                if sn_type=="IIb":
                    where_peak = np.where(np.logical_and(xs>10.,xs<45.))[0]
                else:
                    where_peak = np.where(xs<45.)[0]
                
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
                Ni_Khatami = Ni_K(10**Lp,tp,beta=4/3.)
                ni_t = np.arange(np.argmin([tp-20,0]),np.max(t)+5,0.01)
                print tp,Lp,M_p,Ni_Khatami
            
                Ks = []
                betas = np.arange(0.6,2.1,0.05)
                for beta in betas:
                
                    try:
                        Ks.append(Ni_K(10**Lp,tp,beta=beta,Ni_min=0.005,Ni_max = 1.0)/M_p)
                    except:
                        print "failed at %s"%beta
                        Ks.append(0.0)

                pl.plot(betas,Ks,label=r'$log(L_{peak}(t_{peak}=%2.2f) = %2.2f)$'%(tp,Lp))
        
            except:
                print "No directory %s/"%SN            
                
        except:
            os.chdir(current_dir)
            print "No directory %s/"%SN

    pl.axvline(4/3.,linestyle='--',color='k')
    pl.axhline(1,linestyle='--',color='k')
    for beta_SN,beta_label in [(9/8.,"Ibc"),(0.82,"IIb/pec"),(2.0,'Generic/Mixed'),(0.7,'Central/Helium')]:
        pl.axvline(beta_SN,linestyle='--',linewidth=4,color='k',alpha=0.7)
        pl.annotate(r"$\mathrm{%s}$"%(beta_label),(beta_SN-0.07,1.3),size=20,rotation=90)

    pl.xlabel(r'$\beta$')
    pl.ylabel(r'$M(^{56}\mathrm{Ni})$')
    pl.legend(loc='best',prop={'size':10},ncol=2)
    pl.show()
    os.chdir(current_dir)

    def Ni_K_A(beta,tpeak):
        
        t_ni = 8.8        # niquel decay [days]                                          
        t_co = 111.3      # Cobalt decay [days]                                            
        x = beta*tpeak
        x_ni = x/t_ni
        x_co = x/t_co
        e_ni = 6.395*1e43
        e_co = 1.362*1e43
        corcho1 = e_ni*t_ni*(t_ni*(1-np.exp(-x_ni)) - x*np.exp(-x_ni))
        corcho2 = e_co*t_co*(t_co*(1-np.exp(-x_co)) - x*np.exp(-x_co))
        return x**2. * Q_t(tpeak,M_ni=1.0) / 2. / (corcho1+corcho2)

    tps = np.arange(10,30,1)
    X,Y = np.meshgrid(betas, tps) # grid of point
    Z = Ni_K_A(X,Y)
    im = pl.imshow(Z,cmap=pl.cm.winter) # drawing the function
    # adding the Contour lines with labels
    cset = pl.contour(Z,np.arange(0.5,1.5,0.1),linewidths=2,cmap=pl.cm.Wistia,origin='lower')
    pl.clabel(cset,inline=True,fmt='%1.1f',fontsize=14)
    cb = pl.colorbar(im) # adding the colobar on the right
    cb.set_label(r'$\mathrm{M(Ni_{Khatami})/M(Ni_{Arnett})}$',size=20)
    nx = betas.shape[0]
    no_labels = 5 # how many labels to see on axis x
    step_x = int(nx / (no_labels - 1)) # step between consecutive labels
    x_positions = np.arange(0,nx,step_x) # pixel count at label position
    x_labels = ["%2.2f"%b for b in betas[::step_x]] # labels you want to see
    pl.xticks(x_positions, x_labels)
    
    ny = tps.shape[0]
    no_labels = 5 # how many labels to see on axis x
    step_y = int(ny / (no_labels - 1)) # step between consecutive labels
    y_positions = np.arange(0,ny,step_y) # pixel count at label position
    y_labels = ["%2.2f"%b for b in tps[::step_y]] # labels you want to see
    pl.yticks(y_positions, y_labels)
    pl.xlabel(r"$\beta$",size=25)
    pl.ylabel(r"$t_{peak}$",size=25)
    pl.tight_layout()
    pl.show()
    
if check_t_0:
    
    filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(SN_DATA))]       
    pl.gca().set_color_cycle(filters_colors)
    
    for SN,SN_type,z_SN,E_B_V,t_0,d_L,t_disc,m_disc,t_non,m_non in SN_DATA[['sn','type','sn_z','sn_ebv','t_0','hostlumdist','t_discov','m_discov','t_non_det','m_non_det']]:
        band_string = ("BVRIYJH","BVRIYJH")
        #band_string = ("BVRI","BVRcI")
        print "####### %s ########## \n"%SN
        try: 
            os.chdir("%s/"%SN)
            t_half_1,t_half_2= None,None
            try:
                
                try:
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[0])).T
                    band_string = band_string[0]
                except:
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[1])).T
                    band_string = band_string[1]


                L_disc = 4*pi*(pc_to_cm*d_L*1e6)**2. * 10**(-0.4*(m_disc+21.1))
                L_non = 4*pi*(pc_to_cm*d_L*1e6)**2. * 10**(-0.4*(m_non+21.1))    
                print L_disc,L_non
                xs,ys = t[np.argsort(t)],logL[np.argsort(t)]
                t,logL = xs,ys
                #pl.plot(t,logL,marker='o',color='k')
                #pl.axvline(t_non-t_0)
                #pl.axvline(t_disc-t_0)
                
                where_peak = np.where(np.logical_and(xs>10.,xs<45.))[0]
                xs,ys = xs[where_peak],(10**ys[where_peak])*1e-42
                #k0 = smooth.NonParamRegression(xs, ys, method=npr_methods.SpatialAverage())
                k0 = smooth.NonParamRegression(xs, ys, method = npr_methods.LocalPolynomialKernel(q=4))
                if len(t)>9:
                    k0_full = smooth.NonParamRegression(t,(10**logL)*1e-42, method = npr_methods.LocalPolynomialKernel(q=7))
                else:
                    k0_full = smooth.NonParamRegression(t,(10**logL)*1e-42, method = npr_methods.LocalPolynomialKernel(q=3))
        
                k0.fit()
                k0_full.fit()
                grid = np.arange(np.min(xs),np.max(xs),0.1)
                grid_full = np.arange(np.min(t),np.max(xs),0.1)
                y_grid = np.log10(k0(grid))+42
                y_grid,grid = y_grid[~np.isnan(y_grid)],grid[~np.isnan(y_grid)]
                y_grid_full = np.log10(k0_full(grid_full))+42
                y_grid_full,grid_full = y_grid_full[~np.isnan(y_grid_full)],grid_full[~np.isnan(y_grid_full)]
                where_peak = np.argmax(y_grid)
                tp_0,Lp_0 = grid[where_peak],y_grid[where_peak]
                M_p_0 = M_ni(tp_0,10**Lp_0)
                Ni_Khatami_0 = Ni_K(10**Lp_0,tp_0,beta=4/3.)
                print tp_0,Lp_0,M_p_0,Ni_Khatami_0
                delta_ts = np.linspace(-0.5*(t_disc-t_non),0.5*(t_disc-t_non),10)

                delta_nis = [[],[]]
                
                t = (t+t_0)
                dt = (t_non+t_disc)/2.
                t = t-dt
                Lp = Lp_0
                tp = tp_0+t_0-dt
                M_p_new = M_ni(tp,10**Lp)
                Ni_Khatami_new = Ni_K(10**Lp_0,tp,beta=4/3.)
                print tp,Lp,M_p_new,Ni_Khatami_new
                
                print "old,new explosion epoch : %s,%s"%(t_0,dt)
                print "moving (-%s,%s)"%(0.5*(t_disc-t_non),0.5*(t_disc-t_non))


                for delta_t in delta_ts:
                    
                    try:
                        
                        M_p = M_ni(tp+delta_t,10**Lp)
                        Ni_Khatami = Ni_K(10**Lp,tp+delta_t,beta=4/3.)
                        #print M_p,Ni_Khatami
                        delta_nis[0].append(M_p)
                        delta_nis[1].append(Ni_Khatami)

                    except:
                        
                        delta_nis[0].append(M_p_0)
                        delta_nis[1].append(Ni_Khatami_0)
                        
                delta_nis = np.asarray(delta_nis)
                plot = pl.plot(delta_ts,(delta_nis[0]-M_p_0)/M_p_0,label="%s (%2.2f)"%(SN,tp_0),marker='s',linestyle='--')
                color = plot[0].get_color()
                pl.plot(delta_ts,(delta_nis[1]-Ni_Khatami_0)/(Ni_Khatami_0),color=color,marker='o',linestyle='--')
                pl.axvline(tp_0-tp,color=color)
                
                os.chdir(current_dir)
            except:
                os.chdir(current_dir)
                print "Failed"
                
        except:
            os.chdir(current_dir)
            print "Failed"


            
    pl.axhline(0.0,color='k')
    pl.axvline(0.0,color='k')
    pl.legend(loc='best',prop={'size':10},ncol=2)
    pl.show()
                
    trs = np.arange(7.5,30,2.5)
    filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(trs))]       
    pl.gca().set_color_cycle(filters_colors)
    delta_ts = np.arange(-7,8,1)
    for tr in trs:

        pl.plot(delta_ts,Q_t(tr,1)/Q_t(tr-delta_ts,1),label=r"$t_{rise} \ =\ %2.1f $"%tr)
        
    pl.legend(loc='best',prop={"size":15},ncol=2) 
    pl.xlabel(r"$\Delta t_0$",size=18)
    pl.ylabel(r"$M(^{56}\mathrm{Ni})(t_{rise}-\Delta t_0)/M(^{56}\mathrm{Ni})(t_{rise})$",size=18)
    pl.show()

SN_plot = {'IIb':{"marker":"o","color":"b",'string':'$\mathrm{IIb}$'},'Ib':{"marker":"p","color":"g",'string':'$\mathrm{Ib}$'},'Ic':{"marker":"s","color":"k",'string':'$\mathrm{Ic}$'},'Ibc':{"marker":"x","color":"k",'string':'$\mathrm{Ibc}$'},'Ic_GRB':{"marker":"*","color":"r",'string':'$\mathrm{Ic-GRB}$'},'Ic_BL':{"marker":"s","color":"r",'string':'$\mathrm{Ic_{BL}}$'}}

Ni_0,Ni_K0 = [],[]

fig1 = pl.figure()
ax1 = fig1.add_subplot(111)
fig2 = pl.figure()
ax2 = fig2.add_subplot(111)


if check_type:

    
    for SN,SN_type,z_SN,E_B_V,t_0,d_L,t_disc,m_disc,t_non,m_non in SN_DATA[['sn','type','sn_z','sn_ebv','t_0','hostlumdist','t_discov','m_discov','t_non_det','m_non_det']]:
        band_string = ("BVRIYJH","BVRIYJH")
        #band_string = ("BVRI","BVRcI")
        print "####### %s ########## \n"%SN
        try: 
            os.chdir("%s/"%SN)
            t_half_1,t_half_2= None,None
            try:
                
                try:
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[0])).T
                    band_string = band_string[0]
                except:
                    t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string[1])).T
                    band_string = band_string[1]
                
                if SN_type in ["Ib","Ic","Ibc","Ic_GRB","Ic_BL"]:
                    beta_SN = 9/8.
                else:
                    beta_SN = 0.82

                xs,ys = t[np.argsort(t)],logL[np.argsort(t)]
                t,logL = xs,ys

                if SN!="SN2013cu":
                    where_peak = np.where(np.logical_and(xs>8.,xs<45.))[0]
                else:
                    where_peak = np.where(np.logical_and(xs>0.))[0]

                xs,ys = xs[where_peak],(10**ys[where_peak])*1e-42
                #k0 = smooth.NonParamRegression(xs, ys, method=npr_methods.SpatialAverage())
                k0 = smooth.NonParamRegression(xs, ys, method = npr_methods.LocalPolynomialKernel(q=4))
                
                if len(t)>9:
                    k0_full = smooth.NonParamRegression(t,(10**logL)*1e-42, method = npr_methods.LocalPolynomialKernel(q=7))
                else:
                    k0_full = smooth.NonParamRegression(t,(10**logL)*1e-42, method = npr_methods.LocalPolynomialKernel(q=3))
                
                k0.fit()
                k0_full.fit()
                grid = np.arange(np.min(xs),np.max(xs),0.1)
                grid_full = np.arange(np.min(t),np.max(xs),0.1)
                y_grid = np.log10(k0(grid))+42
                y_grid,grid = y_grid[~np.isnan(y_grid)],grid[~np.isnan(y_grid)]
                y_grid_full = np.log10(k0_full(grid_full))+42
                y_grid_full,grid_full = y_grid_full[~np.isnan(y_grid_full)],grid_full[~np.isnan(y_grid_full)]
                where_peak = np.argmax(y_grid)
                tp_0,Lp_0 = grid[where_peak],y_grid[where_peak]
                M_p_0 = M_ni(tp_0,10**Lp_0)
                Ni_Khatami_0 = Ni_K(10**Lp_0,tp_0,beta=beta_SN)
                Ni_0.append(M_p_0)
                Ni_K0.append(Ni_Khatami_0)
                print SN_type,tp_0,Lp_0,M_p_0,Ni_Khatami_0
                
                ax1.plot(M_p_0,Ni_Khatami_0,marker=SN_plot[SN_type]["marker"],color=SN_plot[SN_type]["color"],markersize=9,alpha=0.9)
                ax2.plot(np.log10(M_p_0),np.log10(M_p_0)-np.log10(Ni_Khatami_0),marker=SN_plot[SN_type]["marker"],color=SN_plot[SN_type]["color"],markersize=9,alpha=0.9)
                os.chdir(current_dir)

                if SN_plot[SN_type].has_key("N"):
                    SN_plot[SN_type]["N"] += 1
                else:
                    SN_plot[SN_type]["N"] = 1
            except:
                os.chdir(current_dir)
                print "Failed"
                
        except:
            os.chdir(current_dir)
            print "Failed"
        
    for type in SN_plot:
        
        if SN_plot[type].has_key("N"):
            ax1.plot(-1,-1,marker=SN_plot[type]["marker"],color=SN_plot[type]["color"],label=r"%s"%SN_plot[type]["string"],linestyle='None',markersize=10)
            ax2.plot(-10,-10,marker=SN_plot[type]["marker"],color=SN_plot[type]["color"],label=r"%s"%SN_plot[type]["string"],linestyle='None',markersize=10)

    ax1.legend(loc='best',prop={'size':18})
    ax2.legend(loc='best',prop={'size':18})
    ax1.set_xlabel(r"$\mathrm{Arnett} \ ^{56}Ni \ [M_\odot]$",size=15)
    ax2.set_xlabel(r"$\mathrm{log Arnett} \ ^{56}Ni \ [M_\odot]$",size=15)
    ax1.set_ylabel(r"$\mathrm{Khatami} \ ^{56}Ni \ [M_\odot]$",size=15)
    ax2.set_ylabel(r"$\mathrm{log Khatami} \ ^{56}Ni \ [M_\odot]$",size=15)
    ax1.plot(np.linspace(0,np.max(Ni_0)+0.01,100),np.linspace(0,np.max(Ni_0)+0.01,100),color='grey')
    #ax2.plot(np.linspace(-1.7,np.log10(np.max(Ni_0)+0.02),100),np.linspace(-1.7,np.log10(np.max(Ni_0)+0.02),100),color='grey')
    ax1.set_xlim(0,np.max(Ni_0)+0.01)
    ax1.set_ylim(0,np.max(Ni_0)+0.01)
    #ax2.set_xlim(-1.7,np.log10(np.max(Ni_0)+0.02))
    #ax2.set_ylim(-1.7,np.log10(np.max(Ni_0)+0.02))
    pl.show()
    
