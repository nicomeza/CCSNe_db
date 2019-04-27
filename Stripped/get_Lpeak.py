from astro_utils import * 
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
sn_labels = ['sn','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0']
sn_formats = ['S15','S20','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_data.dat',dtype={'names':sn_labels,'formats':sn_formats})

compare_nis = False
get_peak = True


M_ni = lambda t_r,L_bol : L_bol/Q_t(t_r,1)

band_string = "BVRI"
current_dir = os.getcwd()

def weird_Ni(t,Mni):
    
    t_ni = 8.8        # niquel decay [days]                                                                                                                            
    t_co = 111.3      # Cobalt decay [days]                                                                                                                            
    Q = Mni*(6.45*np.exp(-t/t_ni)+1.45*np.exp(-t/t_co))*1e43  # [erg/s] . See Nadyozhin 1994   
    Q_dot = -Mni*(6.45*np.exp(-t/t_ni)/t_ni+1.45*np.exp(-t/t_co)/t_co)*1e43
    Weird = Q - Q_dot*t/(2.*t**2. - 1)
    return Weird

if get_peak:
    Lfile = open('Ni56_%s_peak.dat'%band_string,'w')
    Lfile.write("###############%s#######################\n"%band_string)
    Lfile.write("# SN # tpeak \t Lpeak [erg/s] \t log(Lpeak) \t M_ni\n")

    for SN,z_SN,E_B_V,t_0,d_L in SN_DATA[['sn','sn_z','sn_ebv','t_0','hostlumdist']]:
        print "####### %s ########## \n"%SN
        try: 
            os.chdir("%s/"%SN)
            try:
                t,L,logL,Mni = np.genfromtxt('%s_Lbol_%s.dat'%(SN,band_string)).T
                xs,ys = t[np.argsort(t)],logL[np.argsort(t)]
                where_peak = np.where(np.logical_and(xs>10.,xs<40.))[0]
                xs,ys = xs[where_peak],ys[where_peak]
            #k0 = smooth.NonParamRegression(xs, ys, method=npr_methods.SpatialAverage())
                k0 = smooth.NonParamRegression(xs, ys, method = npr_methods.LocalPolynomialKernel(q=3))
                k0.fit()
                grid = np.arange(np.min(xs),np.max(xs),0.1)
                y_grid = k0(grid)
                pl.plot(grid,y_grid,linestyle='--')
                pl.plot(t,logL,marker='o',color='k',linestyle='None',label=SN)

                where_peak = np.argmax(y_grid)
                tp,Lp = grid[where_peak],y_grid[where_peak]
                M_p = M_ni(tp,10**Lp)    
                ni_t = np.arange(np.max([tp-15,0]),np.max(t)+5,0.1)
                pl.plot(ni_t,np.log10(Q_t(ni_t,M_p)),linestyle='--',color='k')
                pl.plot(ni_t,np.log10(weird_Ni(ni_t,M_p)),linestyle='--',color='r')
                pl.xlabel('t-t_0')
                pl.ylabel('log10(L)')
                pl.legend()
                pl.savefig("../%s_%s_Lbol.png"%(SN,band_string))
                pl.show()

                s = "%s \t %s \t %s \t %s \t %s\n"%(SN,tp,10**Lp,Lp,M_p)
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

    nis = np.genfromtxt('Prentice/BVRI_stats.dat',usecols=[0,1,2,3,4,5,6,7,8,9],names=('SN','type','logLp','Lperr1','Lperr2','Mni','Mni_err1','Mni_err2','tp','tp_err'),dtype=('S15','S10','f8','f8','f8','f8','f8','f8','f8','f8'))
    my_nis = np.genfromtxt('Ni56_BVRI_peak.dat',usecols=[0,1,2,3,4],names=('SN','tp','Lp','logLp','Mni'),dtype=('S15','f8','f8','f8','f8'))
    i1,i2 = overlap(my_nis['SN'],nis['SN'])
    
    for item in [('Mni'),('logLp'),('tp')]:
        print item 
        pl.figure()
        pl.scatter(my_nis[i1][item],nis[i2][item],marker='o')
        xmin,xmax = np.min([my_nis[i1][item],nis[i2][item]]),np.max([my_nis[i1][item],nis[i2][item]])
        x = np.arange(xmin-0.3*(xmax-xmin),0.3*(xmax-xmin)+xmax,0.01)
        print xmin,xmax
        pl.xlabel('My %s'%item)
        pl.ylabel('Prentice %s'%item)
        pl.plot(x,x)
        pl.xlim(xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin))
        pl.ylim(xmin-0.1*(xmax-xmin),xmax+0.1*(xmax-xmin))
        
        for sn,mni,mni2 in zip(my_nis[i1]['SN'],my_nis[i1][item],nis[i2][item]):

            pl.annotate(sn,(mni+0.01,mni2+0.01),size=8)
        pl.savefig("%s_comp.png"%item)
    pl.show()
