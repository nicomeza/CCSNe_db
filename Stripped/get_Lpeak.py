from astro_utils import * 
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
from scipy.misc import derivative
from scipy.optimize import newton,brentq

sn_labels = ['sn','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0']
sn_formats = ['S15','S20','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_data.dat',dtype={'names':sn_labels,'formats':sn_formats})

compare_nis = True
get_peak = True
plot = True
show = False

M_ni = lambda t_r,L_bol : L_bol/Q_t(t_r,1)

#band_string = "W1UBVRIJHKs"
band_string = "BVRI"
current_dir = os.getcwd()

def weird_Ni(t,Mni,t_m=1.0):
    
    t_ni = 8.8        # niquel decay [days]                                                                                                                          
    t_co = 111.3      # Cobalt decay [days]                                                           
    Q = Mni*(6.45*np.exp(-t/t_ni)+1.45*np.exp(-t/t_co))*1e43  # [erg/s] . See Nadyozhin 1994   
    Q_dot = -Mni*(6.45*np.exp(-t/t_ni)/t_ni+1.45*np.exp(-t/t_co)/t_co)*1e43
    Weird = Q - t*Q_dot/((t/t_m)**2. - 1)
    return Weird

def weird_Ni_der(t,Mni,t_m=1.0):
    
    t_ni = 8.8        # niquel decay [days]                                                                                          
    t_co = 111.3      # Cobalt decay [days]                                                                                   
    Q = Mni*(6.45*np.exp(-t/t_ni)+1.45*np.exp(-t/t_co))*1e43  # [erg/s] . See Nadyozhin 1994   
    Q_dot = -Mni*(6.45*np.exp(-t/t_ni)/t_ni+1.45*np.exp(-t/t_co)/t_co)*1e43
    x = t/t_m
    Weird = 2*x**2. * Q_dot/(2.*x**2. - 1)
    return Weird
    
def Q_t_dot(t,Mni):
    t_ni = 8.8        # niquel decay [days]                                                                                               
    t_co = 111.3      # Cobalt decay [days]
    return -Mni*(6.45*np.exp(-t/t_ni)/t_ni+1.45*np.exp(-t/t_co)/t_co)*1e43


inflection = True
pl.close()
if get_peak:

    Lfile = open('Ni56_%s_peak.dat'%band_string,'w')
    Lfile.write("###############%s#######################\n"%band_string)
    Lfile.write("# SN # tpeak \t Lpeak [erg/s] \t log(Lpeak) \t M_ni \t t_1/2(<tp) \t t_1/2(>tp) \t FWHM \t t_inflection \t L_dot(inflection)\n")
    for SN,z_SN,E_B_V,t_0,d_L in SN_DATA[['sn','sn_z','sn_ebv','t_0','hostlumdist']]:
        print "####### %s ########## \n"%SN
        try: 
            os.chdir("%s/"%SN)
            t_half_1,t_half_2= None,None
            try:
                
                t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string)).T
                xs,ys = t[np.argsort(t)],logL[np.argsort(t)]
                t,logL = xs,ys
                where_peak = np.where(np.logical_and(xs>10.,xs<45.))[0]
                xs,ys = xs[where_peak],(10**ys[where_peak])*1e-42
                #k0 = smooth.NonParamRegression(xs, ys, method=npr_methods.SpatialAverage())
                k0 = smooth.NonParamRegression(xs, ys, method = npr_methods.LocalPolynomialKernel(q=4))
                k0_full = smooth.NonParamRegression(t,(10**logL)*1e-42, method = npr_methods.LocalPolynomialKernel(q=4))
                k0.fit()
                k0_full.fit()
                grid = np.arange(np.min(xs),np.max(xs),0.1)
                grid_full = np.arange(np.min(t),np.max(xs),0.1)
                y_grid = np.log10(k0(grid))+42
                y_grid = y_grid[~np.isnan(y_grid)]
                y_grid_full = np.log10(k0_full(grid_full))+42
                y_grid_full = y_grid_full[~np.isnan(y_grid_full)]
                where_peak = np.argmax(y_grid)
                tp,Lp = grid[where_peak],y_grid[where_peak]
                M_p = M_ni(tp,10**Lp)    
                ni_t = np.arange(np.argmin([tp-20,0]),np.max(t)+5,0.01)
                print tp,Lp,M_p
                
                if plot:

                    pl.figure()
                    pl.plot(grid_full,y_grid_full,linestyle='--',color='k')
                    pl.plot(grid,y_grid,color='g')
                    pl.plot(t,logL,marker='o',color='k',linestyle='None',label=SN,alpha=0.6)
                    pl.plot(ni_t,np.log10(Q_t(ni_t,M_p)),linestyle='-.',color='k')
                    pl.axhline(Lp-np.log10(2),linestyle='--',color='k')
                    pl.annotate("Half max",(np.max(t)-5,Lp-np.log10(2)+0.05),size=8)
                    
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
                    der2 = derivative(k0_full, grid_full[before_p_full], dx=0.5, n=2, order=order)
                    inf_t = grid_full[before_p_full][np.argmin(np.abs(der2))]
                    der_at_inf =  1e42*derivative(k0, inf_t, dx=0.5, n=2, order=order)
                    print (der_at_inf - Q_t_dot(inf_t,M_p))/Q_t_dot(inf_t,M_p)
                    
                    if plot:
                        for t_m in [5,10,15,20,25,30]:
                            print t_m,(weird_Ni_der(inf_t,M_p,t_m=t_m) - der_at_inf)/der_at_inf
                        pl.axvline(inf_t,color='k',linestyle='--')
                        pl.annotate("Inflection",(inf_t-2.,np.min(logL)+0.1),size=8,rotation='vertical')
                if len(after_p):
                    print "After peak :\n"
                    try:
                        
                        f = lambda x : k0_full(x) - 10**(Lp-42)/2.
                        zero = brentq(f,tp+2,np.max(t))     
                        t_half_2 = zero - tp
                        print "t_1/2 = %s"%t_half_2
                    except:
                        print "No t_1/2"

                    der2 = derivative(k0, grid[after_p], dx=0.5, n=2, order=order)
                    inf_t2 = grid[after_p][np.argmin(np.abs(der2))]
                    der_at_inf_2 =  1e42 * derivative(k0, inf_t2, dx=0.5, n=2, order=order)
                    print (der_at_inf_2 - Q_t_dot(inf_t2,M_p))/Q_t_dot(inf_t2,M_p)
                    for t_m in [5,10,15,20,25,30]:
                        print t_m,(weird_Ni_der(inf_t2,M_p,t_m=t_m) - der_at_inf_2)/der_at_inf_2                    
                    if plot:
                        pl.axvline(inf_t2,color='k',linestyle='--')
                        pl.annotate("Inflection",(inf_t2-2.,np.min(logL)+0.1),size=8,rotation='vertical')
                if plot:
                    if inflection:
                        for t_m,color in zip([5,10,20,30],['g','r','b','orange']):
                            pl.plot(ni_t,np.log10(weird_Ni(ni_t,M_p,t_m=t_m)),linestyle='--',color=color,label='t_m = %s'%t_m)                    
                    
                try:
                    FWHM = t_half_2+t_half_1
                    print "FWHM: %s"%FWHM
                except:
                    FWHM = 99
                if plot:
                    pl.xlabel('t-t_0')
                    pl.ylabel('log10(L)')
                    pl.legend()
                    pl.ylim(np.min(logL)-0.1,Lp+0.2)
                    pl.savefig("../%s_%s_Lbol.png"%(SN,band_string))
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
                        pl.close()
                    
                s = "%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n"%(SN,tp,10**Lp,Lp,M_p,t_half_1,t_half_2,FWHM,inf_t2,der_at_inf_2)
                print s 
                Lfile.write(s)
                os.chdir(current_dir)
            except:
                print "No Lbol file"
                os.chdir(current_dir)
        except:
            print "No directory %s/"%SN
        
        
    Lfile.close()



if compare_nis:

    pl.close()
    nis = np.genfromtxt('Prentice/BVRI_stats.dat',\
                            names=('SN','type','logLp','logLp_err1','logLp_err2','Mni','Mni_err1','Mni_err2','tp','tp_err','trise','trise_err','tdecay','tdecay_err','FWHM','FWHM_err'),\
                            dtype=('S15','S10','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))
    my_nis = np.genfromtxt('Ni56_BVRI_peak.dat',\
                               names=('SN','tp','Lp','logLp','Mni','trise','tdecay','FWHM','t_inf','L_dot_inf'),\
                               dtype=('S15','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

    i1,i2 = overlap(my_nis['SN'],nis['SN'])
    
    for item in [('Mni'),('logLp'),('tp'),('trise'),('tdecay'),('FWHM')]:
        print item 
        pl.figure()
        no_99 = np.logical_and(nis[i2][item] != 99.,my_nis[i1][item] != 99.)
        no_nan = np.logical_and(~np.isnan(my_nis[i1][item]),no_99)
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
        pl.savefig("%s_comp.png"%item)
    if show:
        pl.show()

    for item1,item2 in [('FWHM','L_dot_inf')]:
        no_99 = np.logical_and(my_nis[item1] != 99.,my_nis[item2] != 99.)
        no_nan = np.logical_and(~np.isnan(my_nis[item1]),~np.isnan(my_nis[item2]))
        valid = np.where(np.logical_and(no_99,no_nan))[0]
        pl.plot(my_nis[item1][valid],my_nis[item2][valid],'o')
    pl.show()
        
