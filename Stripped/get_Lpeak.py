from astro_utils import * 
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
sn_labels = ['sn','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0']
sn_formats = ['S15','S20','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_data.dat',dtype={'names':sn_labels,'formats':sn_formats})

compare_nis = True
get_peak = False 

if compare_nis:

    nis = np.genfromtxt('Prentice/BVRI_stats.dat',usecols=[0,1,2,3,4,5,6,7,8,9],names=('SN','type','Lp','Lperr1','Lperr2','Mni','Mni_err1','Mni_err2','tp','tp_err'),dtype=('S15','S10','f8','f8','f8','f8','f8','f8','f8','f8'))
    my_nis = np.genfromtxt('Ni56_BVRI_peak.dat',usecols=[0,1,2,3,4],names=('SN','tp','Lp','logLp','Mni'),dtype=('S15','f8','f8','f8','f8'))
    i1,i2 = overlap(my_nis['SN'],nis['SN'])
    pl.plot(my_nis[i1]['Mni'],nis[i2]['Mni'],marker='o',linestyle='None')
    x = np.arange(0,0.5,0.01)
    pl.xlabel('My Ni56')
    pl.ylabel('Prentice Ni56')
    pl.plot(x,x)

    pl.figure()
    pl.plot(np.log10(my_nis[i1]['Lp']),nis[i2]['Lp'],marker='o',linestyle='None')
    x = np.arange(41.1,43.5,0.1)
    pl.xlabel('My Lp')
    pl.ylabel('Prentice Lp')
    pl.plot(x,x)

    pl.figure()
    pl.plot((my_nis[i1]['tp']),nis[i2]['tp'],marker='o',linestyle='None')
    x = np.arange(8,40,1)
    pl.xlabel('My tp')
    pl.ylabel('Prentice tp')
    pl.plot(x,x)

    pl.show()

M_ni = lambda t_r,L_bol : L_bol/Q_t(t_r,1)

band_string = "BVRI"
current_dir = os.getcwd()

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
                pl.plot(xs,ys,marker='o',color='k')
                pl.show()
                where_peak = np.argmax(y_grid)
                tp,Lp = grid[where_peak],y_grid[where_peak]
                M_p = M_ni(tp,10**Lp)    
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
