import numpy as np
from scipy import stats
import matplotlib.pyplot as pl

sn_labels = ['sn','type','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_test_total.dat',dtype={'names':sn_labels,'formats':sn_formats})


SN_plot = {'IIb':{"marker":"o","color":"b",'string':'$\mathrm{IIb}$'},\
               'Ib':{"marker":"p","color":"g",'string':'$\mathrm{Ib}$'},\
               'Ic':{"marker":"s","color":"k",'string':'$\mathrm{Ic}$'},\
               'Ibc':{"marker":"x","color":"k",'string':'$\mathrm{Ibc}$'},\
               'Ic_GRB':{"marker":"*","color":"r",'string':'$\mathrm{Ic-GRB}$'},\
               'Ic_BL':{"marker":"s","color":"r",'string':'$\mathrm{Ic_{BL}}$'}}


colormap = pl.cm.spectral

#my_file = "Ni56_BVRIYJH_local_peak.dat"
#save_name = "BVRIYJH_local"

my_file = "Ni56_BVRIYJH_total_peak.dat"
save_name = "BVRIYJH_total"

my_nis = np.genfromtxt(my_file,\
                           names=('SN','type','tp','Lp','logLp','Mni_peak','Mni_Khatami','Mni_tail','trise','tdecay','FWHM','t_inf','L_dot_inf'),\
                           dtype=('S15','S10','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

nis_II = [l.split() for l in open('SNII_Ni56.dat').readlines()[1:]]
nis_II = [(sn,np.median(np.asarray(n.split(','),dtype='f8'))) for sn,n in nis_II]


Ni_types = ['Mni_tail','Mni_Khatami','Mni_peak']
sorted_data_II = np.sort([x[1] for x in nis_II])


fig2 = pl.figure()
ax2 = fig2.add_subplot(111)
ax2.step(np.concatenate([sorted_data_II,[0]]),np.arange(0.0,sorted_data_II.size+1)/sorted_data_II.size,\
            color='r',label="II (%s)"%len(nis_II),linewidth=2)

for ni_type,sn_type in zip(Ni_types,['Ib','IIb','Ic']):

    print "#############%s##################"%ni_type
    KS_test,AD_test = [],[]
    fig = pl.figure()
    ax = fig.add_subplot(111)
    where_type = np.where(np.logical_and(np.logical_and(my_nis['type']!="Ic_BL",my_nis['type']!="Ic_GRB"),~np.isnan(my_nis['%s'%ni_type])))[0]
    if len(where_type)>1:
        #print "Number of %s : %s"%(sn_type,len(where_type))                                                                                                                                                                                                                
        sorted_data = np.sort(my_nis['%s'%ni_type][where_type])
        ax.step(np.concatenate([sorted_data,[0]]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,\
                    color=SN_plot[sn_type]['color'],label="IIb+Ibc (%s)"%len(where_type),linewidth=2)
        ax2.step(np.concatenate([sorted_data,[0]]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,\
                    color=SN_plot[sn_type]['color'],label="IIb+Ibc : %s (%s)"%(ni_type.split('_')[1],len(where_type)),linewidth=2)

        mask = [True]*len(where_type)
        for index in np.arange(0,len(where_type)):
            mask[index] = False
            ks = stats.ks_2samp(my_nis['%s'%ni_type][where_type[mask]],sorted_data_II)
            ad = stats.anderson_ksamp([my_nis['%s'%ni_type][where_type][mask],sorted_data_II])
            if ~np.isnan(ks[1]) and ~np.isnan(ks[0]):
                KS_test.append(ks)
            else:
                print "NaN value at %s"%my_nis['SN'][where_type[index]]
            if ~np.isnan(ad[2]) and ~np.isnan(ad[0]):
                AD_test.append(ad)
            else:
                print "NaN value at %s"%my_nis['SN'][where_type[index]]
            mask[index] = True
            
    print "### Kolmogorov-Smirnov ###"
    D,KS = [ks[0] for ks in KS_test],[ks[1] for ks in KS_test]
    KS_mean,KS_std =  np.mean(KS),np.std(KS)
    D_mean,D_std =  np.mean(D),np.std(D)
    print "Full data : %s,%s"%stats.ks_2samp(my_nis['%s'%ni_type][where_type],sorted_data_II)
    D_interval = stats.t.interval(0.99,len(KS)-1,loc=D_mean,scale=D_std)
    KS_interval = stats.t.interval(0.99,len(KS)-1,loc=KS_mean,scale=KS_std)
    print "KS \t D "
    if KS_mean>1e-3:
        print " %3.3f  (%3.3f,%3.3f) & %3.3f (%3.3f, %3.3f) \\"%(KS_mean,(KS_interval)[0],(KS_interval)[1],D_mean,D_interval[0],D_interval[1])
    else:
        print " %3.3E  (%3.3E,%3.3E) & %3.3f (%3.3f,%3.3f) \\"%(KS_mean,(KS_interval)[0],(KS_interval)[1],D_mean,D_interval[0],D_interval[1])
    print "### Anderson-Darling ###"
    S,AD = [ad[0] for ad in AD_test],[ad[2] for ad in AD_test]
    AD_mean,AD_std =  np.mean(AD),np.std(AD)
    S_mean,S_std =  np.mean(S),np.std(S)
    print "Full data : %s,%s,%s"%stats.anderson_ksamp([my_nis['%s'%ni_type][where_type],sorted_data_II])
    S_interval = stats.t.interval(0.99,len(AD)-1,loc=S_mean,scale=S_std)
    AD_interval = stats.t.interval(0.99,len(AD)-1,loc=AD_mean,scale=AD_std)
    print "AD \t S "
    if AD_mean>1e-3:
        print " %3.3f  (%3.3f,%3.3f) & %3.3f (%3.3f, %3.3f) \\"%(AD_mean,(AD_interval)[0],(AD_interval)[1],S_mean,S_interval[0],S_interval[1])
    else:
        print " %3.3E  (%3.3E,%3.3E) & %3.3f (%3.3f,%3.3f) \\"%(AD_mean,(AD_interval)[0],(AD_interval)[1],S_mean,S_interval[0],S_interval[1])

    ax.step(np.concatenate([sorted_data_II,[0]]),np.arange(0.0,sorted_data_II.size+1)/sorted_data_II.size,\
                color='r',label="II (%s)"%len(nis_II),linewidth=2)

    ax.set_title(r"$^{56}Ni \ %s$"%ni_type.split('_')[1],size=15)
    ax.legend(loc='best')
    ax.set_xlabel(r"$\mathrm{Nickel \ Mass} \ M_\odot$",size=15)
    fig.savefig("%s_%s_hist_SESN.png"%(save_name,ni_type))


ax2.legend(loc='best')
ax2.set_xlabel(r"$\mathrm{Nickel \ Mass} \ M_\odot$",size=15)
fig2.savefig("%s_all_hist_SESN.png"%save_name)
pl.show()

