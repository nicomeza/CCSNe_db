from astro_utils import *
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
from scipy.misc import derivative
from scipy.optimize import newton,brentq
from astroML.plotting import hist
from scipy import stats

sn_labels = ['sn','type','host','hostredshift','hostlumdist','sn_ebv','sn_host','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_test_total.dat',dtype={'names':sn_labels,'formats':sn_formats})


SN_plot = {'IIb':{"marker":"o","color":"b",'string':'$\mathrm{IIb}$'},'Ib':{"marker":"p","color":"g",'string':'$\mathrm{Ib}$'},'Ic':{"marker":"s","color":"k",'string'\
:'$\mathrm{Ic}$'},'Ibc':{"marker":"x","color":"k",'string':'$\mathrm{Ibc}$'},'Ic_GRB':{"marker":"*","color":"r",'string':'$\mathrm{Ic-GRB}$'},\
               'Ic_BL':{"marker":"s","color":"r",'string':'$\mathrm{Ic_{BL}}$'}}


colormap = pl.cm.spectral

my_file = "Ni56_BVRIYJH_local_peak.dat"
save_name = "BVRIYJH_local"

my_file = "Ni56_BVRIYJH_total_peak.dat"
save_name = "BVRIYJH_total"

my_file = "Ni56_BVRIYJH_shift_peak.dat"
save_name = "BVRIYJH_shift"

my_nis = np.genfromtxt(my_file,\
                           names=('SN','type','tp','Lp','logLp','Mni_peak','Mni_Khatami','Mni_tail','trise','tdecay','FWHM','t_inf','L_dot_inf'),\
                           dtype=('S15','S10','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

nis_II = [l.split() for l in open('SNII_Ni56.dat').readlines()[1:]]
nis_II = [(sn,np.median(np.asarray(n.split(','),dtype='f8'))) for sn,n in nis_II]

where_tail = np.where(~np.isnan(my_nis['Mni_tail']))[0]
where_rise = np.where(~np.isnan(my_nis['tp']))[0]
nis_rise = my_nis[where_rise]
nis_compare = my_nis[where_tail]

one_one = False
sample_dist = False

############## ONE TO ONE PLOTS #######################33
if one_one:
    fig = pl.figure()
    ax = fig.add_subplot(111)

    nis_line = np.arange(0.0,0.5,0.1)

    for sn_type in SN_plot:
        where_type = np.where(nis_compare['type']==sn_type)[0]
        if len(where_type)>0:
            ax.plot(nis_compare['Mni_peak'][where_type],nis_compare['Mni_tail'][where_type],\
                        marker=SN_plot[sn_type]["marker"],color=SN_plot[sn_type]["color"],linestyle='None',label=r"%s"%SN_plot[sn_type]["string"],markersize=10,alpha=0.8)

            ax.plot(nis_line,nis_line,color='k',linestyle='--')
    ax.set_xlabel(r"$\mathrm{Arnett} \ ^{56}Ni \ [M_\odot]$",size=22)
    ax.set_ylabel(r"$\mathrm{Tail} \ ^{56}Ni \ [M_\odot]$",size=22)
    ax.set_xlim(0,np.max(nis_compare['Mni_peak'])+0.01)
    ax.set_ylim(0,np.max(nis_compare['Mni_peak'])+0.01)
    pl.legend(loc='best',prop={'size':17},ncol=2)
    pl.tick_params(labelsize=17)
#pl.subplots_adjust(left=0.12)                                                                                                                                     
#pl.subplots_adjust(right=0.75)                                                                                                                                    
#pl.subplots_adjust(top=0.95)                                                                                                                                      
    pl.tight_layout()
    pl.setp(pl.legend().get_texts(), fontsize='16')
    pl.savefig("Arnett_Tail_%s.png"%save_name,dpi=300)
    pl.show()


    fig = pl.figure()
    ax = fig.add_subplot(111)

    nis_line = np.arange(0.0,0.5,0.1)

    for sn_type in SN_plot:
        where_type = np.where(my_nis['type']==sn_type)[0]
        if len(where_type)>0:
            ax.plot(my_nis['Mni_peak'][where_type],my_nis['Mni_Khatami'][where_type],\
                        marker=SN_plot[sn_type]["marker"],color=SN_plot[sn_type]["color"],linestyle='None',label=r"%s"%SN_plot[sn_type]["string"],markersize=10,alpha=0.8)

    ax.plot(nis_line,nis_line,color='k',linestyle='--')
    ax.set_xlabel(r"$\mathrm{Arnett} \ ^{56}Ni \ [M_\odot]$",size=22)
    ax.set_ylabel(r"$\mathrm{Khatami} \ ^{56}Ni \ [M_\odot]$",size=22)
    ax.set_xlim(0,np.max(my_nis['Mni_peak'])+0.01)
    ax.set_ylim(0,np.max(my_nis['Mni_peak'])+0.01)
    pl.tick_params(labelsize=17)
#pl.subplots_adjust(left=0.12)                                                                                                                                     
#pl.subplots_adjust(right=0.75)                                                                                                                                    
#pl.subplots_adjust(top=0.95)                                                                                                                                      
    pl.tight_layout()
    pl.setp(pl.legend().get_texts(), fontsize='16')
    pl.legend(loc='upper left',prop={'size':17},ncol=2)
    pl.savefig("Arnett_Khatami_%s.png"%(save_name),dpi=300)
    pl.show()


    pl.hist(nis_compare['Mni_peak'],label='Peak',alpha=0.5)
    pl.hist(nis_compare['Mni_tail'],label='Tail',alpha=0.5)
    pl.legend()
    pl.show()

    pl.hist(nis_compare['Mni_Khatami'],label='Khatami',alpha=0.5)
    pl.hist(nis_compare['Mni_tail'],label='Tail',alpha=0.5)
    pl.legend()
    pl.show()

######################### SAMPLE DISTRIBUTIONS #########################3



## redshift distributions ##
if sample_dist:
    for sn_type in SN_plot:
    
        where_type = np.where(SN_DATA['type']==sn_type)[0]
        if len(where_type)>0:
            print sn_type, len(where_type)
            if len(where_type)>3:
                hist(SN_DATA[where_type]['hostredshift'],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='knuth')
            else:
                hist(SN_DATA[where_type]['hostredshift'],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='blocks')
            pl.axvline(np.median(SN_DATA[where_type]['hostredshift']),linestyle='--',color=SN_plot[sn_type]['color'],linewidth=5)
    
    pl.xlabel(r'$\mathrm{Redshift}$',size=22)
    pl.legend(loc='best',prop={'size':17},ncol=2)
    pl.tick_params(labelsize=17)
#pl.subplots_adjust(left=0.12)                                                                                                                                     
#pl.subplots_adjust(right=0.75)                                                                                                                                    
#pl.subplots_adjust(top=0.95)                                                                                                                                      
    pl.tight_layout()
    pl.setp(pl.legend().get_texts(), fontsize='16')
    pl.savefig("z_hist_subtypes.png",dpi=300)
    pl.show()

## reddening distributions ##

    for sn_type in SN_plot:
    
        where_type = np.where(SN_DATA['type']==sn_type)[0]
        if len(where_type)>0:
            print sn_type, len(where_type)
            if len(where_type)>3:
                hist(SN_DATA[where_type]['sn_ebv'],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='knuth')
            else:
                hist(SN_DATA[where_type]['sn_ebv'],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='blocks')
            pl.axvline(np.median(SN_DATA[where_type]['sn_ebv']),linestyle='--',color=SN_plot[sn_type]['color'],linewidth=5)
    
    pl.legend(loc='best',prop={'size':17},ncol=2)
    pl.tick_params(labelsize=17)
#pl.subplots_adjust(left=0.12)                                                                                                                                     
#pl.subplots_adjust(right=0.75)                                                                                                                                    
    pl.subplots_adjust(bottom=0.05)                                                                                                                                      
    pl.tight_layout()
    pl.setp(pl.legend().get_texts(), fontsize='16')
    pl.xlabel(r'$E(B-V)$',size=22)
    pl.savefig("MW_EBV_hist_subtypes.png",dpi=300)
    pl.show()

    for sn_type in SN_plot:
    
        where_type = np.where(np.logical_and(SN_DATA['type']==sn_type,SN_DATA['sn_host']<10))[0]
        if len(where_type)>0:
            print sn_type, len(where_type)
            if len(where_type)>3:
                hist(SN_DATA[where_type]['sn_host'],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='knuth')
            else:
                hist(SN_DATA[where_type]['sn_host'],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='blocks')
            pl.axvline(np.median(SN_DATA[where_type]['sn_host']),linestyle='--',color=SN_plot[sn_type]['color'],linewidth=5)
    
    pl.legend(loc='best',prop={'size':17},ncol=2)
    pl.tick_params(labelsize=17)
#pl.subplots_adjust(left=0.12)                                                                                                                                     
#pl.subplots_adjust(right=0.75)                                                                                                                                    
#pl.subplots_adjust(top=0.95)                                                                                                                                      
    pl.tight_layout()
    pl.setp(pl.legend().get_texts(), fontsize='16')
    pl.xlabel(r'$E(B-V)$',size=22)
    pl.subplots_adjust(bottom=0.07)
    pl.savefig("host_EBV_hist_subtypes.png",dpi=300)
    pl.show()


## deltat distributions ##

    for sn_type in SN_plot:
    
        where_type = SN_DATA['type']==sn_type
        where_t_0 = np.logical_and(SN_DATA['t_non_det']>99.,SN_DATA['t_discov']>99.)
        where_type = np.where(np.logical_and(where_t_0,where_type))[0]
        if len(where_type)>0:
            print sn_type, len(where_type)
            if len(where_type)>3:
                hist(SN_DATA[where_type]['t_non_det']-SN_DATA[where_type]['t_discov'],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='knuth')
            else:
                hist(SN_DATA[where_type]['t_non_det']-SN_DATA[where_type]['t_discov'],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='blocks')
            pl.axvline(np.median(SN_DATA[where_type]['t_non_det']-SN_DATA[where_type]['t_discov']),linestyle='--',color=SN_plot[sn_type]['color'],linewidth=5)
        else:

            print "No %s"%sn_type

    pl.xlabel(r'$\Delta t$',size=22)
    pl.legend(loc='best',prop={'size':17},ncol=2)
    pl.tick_params(labelsize=17)
#pl.subplots_adjust(left=0.12)                                                                                                                                     
#pl.subplots_adjust(right=0.75)                                                                                                                                    
#pl.subplots_adjust(top=0.95)                                                                                                                                      
    pl.tight_layout()
    pl.setp(pl.legend().get_texts(), fontsize='16')

    pl.savefig("dt_hist_subtypes.png",dpi=300)
    pl.show()



## Rise time distributions ##
    for sn_type in SN_plot:
    
        where_type = np.where(nis_rise['type']==sn_type)[0]
        label_string = ' '.join(sn_type.split('_'))
        if len(where_type)>0:
            print sn_type, len(where_type)
            if len(where_type)>3:
                hist(nis_rise['tp'][where_type],label="%s (%s)"%(label_string,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='knuth')
            else:
            
                hist(nis_rise['tp'][where_type],label="%s (%s)"%(label_string,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='blocks')
            pl.axvline(np.median(nis_rise['tp'][where_type]),linestyle='--',color=SN_plot[sn_type]['color'],linewidth=5)
    
    pl.tick_params(labelsize=17)
    pl.subplots_adjust(left=0.05)                                                                                                                                     
    pl.subplots_adjust(right=0.95)                                                                                                                                    
    pl.subplots_adjust(top=0.98)                                                                                                                                      
    pl.subplots_adjust(bottom=0.12)                                                                                         
    pl.setp(pl.legend().get_texts(), fontsize='16')
    pl.legend(loc='best',prop={'size':17},ncol=2)
    pl.xlabel(r'$\mathrm{Rise \ time \ [days]}$',size=22)
    pl.savefig("rise_hist_subtypes.png",dpi=300)
    pl.show()


dt_ = SN_DATA['t_non_det']-SN_DATA['t_discov']
dt_where = np.where(np.logical_and(SN_DATA['t_discov']>99.,SN_DATA['t_non_det']>99))
dmag_ = SN_DATA['m_non_det']-SN_DATA['m_discov']
dmag_where = np.where(np.logical_and(SN_DATA['m_discov']<99,SN_DATA['m_non_det']<99))

F_discov =4*pi * (SN_DATA['hostlumdist'][dmag_where]*1e6*pc_to_cm)**2. *  10** (-0.4 * (SN_DATA['m_discov'][dmag_where] + 0.21+48.6))  # ASSUMING ALL MAGS ARE FROM R band.
F_non_det =4*pi * (SN_DATA['hostlumdist'][dmag_where]*1e6*pc_to_cm)**2. *  10** (-0.4 * (SN_DATA['m_non_det'][dmag_where] + 0.21+48.6))
F_model = 10** (-0.4 * (-13. + 0.21+48.6)) * 4*pi * (10*pc_to_cm)**2.  # Minimum R mag from model SN5p11Ax1
normed_F_non_det = (F_non_det - F_model)/F_non_det
normed_F_non_det = (F_discov - F_model)/F_discov
abs_mags_discov = SN_DATA['m_discov'][dmag_where] - 5*np.log10(SN_DATA['hostlumdist'][dmag_where]) - 25 # Absolute mags 
abs_mags_non_det = SN_DATA['m_non_det'][dmag_where] - 5*np.log10(SN_DATA['hostlumdist'][dmag_where]) - 25 # Absolute mags 


## Nickel distributions ####

Ni_types = ['Mni_tail','Mni_Khatami','Mni_peak'] 

for ni_type in Ni_types:

    fig = pl.figure()
    ax = fig.add_subplot(111)
    print ni_type
    
    for sn_type in SN_plot:
        label_string = ' '.join(sn_type.split('_'))
        print label_string
        #if sn_type not in ['Ic_BL','Ic_GRB','Ibc','Ib']:
        where_type = np.where(np.logical_and(my_nis['type']==sn_type,~np.isnan(my_nis['%s'%ni_type])))[0]
        if len(where_type)>1:
            print stats.ks_2samp(my_nis['%s'%ni_type][where_type],[x[1] for x in nis_II])
            sorted_data = np.sort(my_nis['%s'%ni_type][where_type]) 
            ax.step(np.concatenate([[0],sorted_data]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,color=SN_plot[sn_type]['color'],label="%s (%s)"%(sn_type,len(where_type)))
            #ax.hist(my_nis['%s'%ni_type][where_type],color=SN_plot[sn_type]['color'],label=sn_type,alpha=0.1,ls='dashed',cumulative=True,normed=True,fill=None)

    ax.set_title(r"$^{56}Ni \ %s$"%ni_type.split('_')[1],size=20)
    pl.tick_params(labelsize=17)
    #pl.subplots_adjust(left=0.12)                                                                                                                                 
    #pl.subplots_adjust(right=0.75)                                                                           
    #pl.subplots_adjust(top=0.95)                                                                                                                                     
    pl.setp(ax.legend().get_texts(), fontsize='16')
    pl.legend(loc='best',prop={'size':17},ncol=2)
    pl.savefig("%s_%s_hist_subtypes.png"%(save_name,ni_type),dpi=300)
    pl.show()


sorted_data_II = np.sort([x[1] for x in nis_II])     
    
for ni_type in Ni_types:

    fig = pl.figure()
    ax = fig.add_subplot(111)
   
    sn_type = "IIb"
    where_type = np.where(np.logical_and(my_nis['type']==sn_type,~np.isnan(my_nis['%s'%ni_type])))[0]
    if len(where_type)>1:
        print "Number of %s : %s"%(sn_type,len(where_type))
        sorted_data = np.sort(my_nis['%s'%ni_type][where_type]) 
        ax.step(np.concatenate([sorted_data,[0]]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,\
                    color=SN_plot[sn_type]['color'],label="IIb (%s)"%(len(where_type)),linestyle='--')
    
    sn_type = "Ic"
    where_type = np.where(np.logical_and(np.logical_and(my_nis['type']!="IIb",my_nis['type']!="Ic_GRB"),~np.isnan(my_nis['%s'%ni_type])))[0]
    if len(where_type)>1:
        print "Number of %s : %s"%(sn_type,len(where_type))
        sorted_data = np.sort(my_nis['%s'%ni_type][where_type]) 
        ax.step(np.concatenate([sorted_data,[0]]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,\
                    color=SN_plot[sn_type]['color'],label="Ibc (%s)"%len(where_type),linestyle='--')
        
    sn_type = "SE-CCSNe"
    where_type = np.where(np.logical_and(np.logical_and(my_nis['type']!="Ic_BL",my_nis['type']!="Ic_GRB"),~np.isnan(my_nis['%s'%ni_type])))[0]
    if len(where_type)>1:
        print "Number of %s : %s"%(sn_type,len(where_type))
        sorted_data = np.sort(my_nis['%s'%ni_type][where_type]) 
        ax.step(np.concatenate([[0],sorted_data]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,\
                    color='magenta',label="IIb/Ibc (%s)"%len(where_type),linewidth=3)
        


    ax.step(np.concatenate([[0],sorted_data_II]),np.arange(0.0,sorted_data_II.size+1)/sorted_data_II.size,\
                color='r',label="II (%s)"%len(nis_II),linewidth=3)
    
    pl.title(r"$^{56}Ni \ %s$"%ni_type.split('_')[1],size=20)
    pl.legend(loc='best',prop={'size':17},ncol=2)
    pl.tick_params(labelsize=17)
#pl.subplots_adjust(left=0.12)                                                                                                                                     
#pl.subplots_adjust(right=0.75)                                                                                                                                    
#pl.subplots_adjust(top=0.95)                                                                                                                                      
    pl.setp(pl.legend().get_texts(), fontsize='16')
    ax.set_xlabel(r"$\mathrm{Nickel \ Mass} \ M_\odot$",size=22)
    pl.savefig("%s_%s_hist_IIb-Ibc.png"%(save_name,ni_type),dpi=300)
    pl.show()




fig2 = pl.figure()
ax2 = fig2.add_subplot(111)   
ax2.step(np.concatenate([sorted_data_II,[0]]),np.arange(0.0,sorted_data_II.size+1)/sorted_data_II.size,\
            color='r',label="II (%s)"%len(nis_II),linewidth=3.5)
    
for ni_type,sn_type,ni_label in zip(Ni_types,['Ib','IIb','Ic'],['Tail','K&K','Arnett\'s']):
    
    print ni_type
    KS_test = []
    fig = pl.figure()
    ax = fig.add_subplot(111)
    where_type = np.where(np.logical_and(np.logical_and(my_nis['type']!="Ic_BL",my_nis['type']!="Ic_GRB"),~np.isnan(my_nis['%s'%ni_type])))[0]
    if len(where_type)>1:
        #print "Number of %s : %s"%(sn_type,len(where_type))
        sorted_data = np.sort(my_nis['%s'%ni_type][where_type]) 
        ax.step(np.concatenate([sorted_data,[0]]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,\
                    color=SN_plot[sn_type]['color'],label="IIb+Ibc (%s)"%len(where_type),linewidth=2)
        ax2.step(np.concatenate([sorted_data,[0]]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,\
                    color=SN_plot[sn_type]['color'],label="IIb+Ibc : %s (%s)"%(ni_label,len(where_type)),linewidth=3.5)
        
        mask = [True]*len(where_type)
        for index in np.arange(0,len(where_type)):
            mask[index] = False
            ks = stats.ks_2samp(my_nis['%s'%ni_type][where_type[mask]],sorted_data_II)
            if ~np.isnan(ks[1]):
                KS_test.append(stats.ks_2samp(my_nis['%s'%ni_type][where_type[mask]],sorted_data_II))
            else:
                print "NaN value at %s"%my_nis['SN'][where_type[index]]
            mask[index] = True
        
    D,KS = [ks[0] for ks in KS_test],[ks[1] for ks in KS_test]   
    KS_mean,KS_std =  np.mean(KS),np.std(KS)
    D_mean,D_std =  np.mean(D),np.std(D)
    print stats.ks_2samp(my_nis['%s'%ni_type][where_type],sorted_data_II)
    D_interval = stats.t.interval(0.99,len(KS)-1,loc=D_mean,scale=D_std)
    KS_interval = stats.t.interval(0.99,len(KS)-1,loc=KS_mean,scale=KS_std)
    print "KS \t D "
    if KS_mean>1e-3:
        print " %3.3f  (%3.3f,%3.3f) & %3.3f (%3.3f, %3.3f) \\"%(KS_mean,(KS_interval)[0],(KS_interval)[1],D_mean,D_interval[0],D_interval[1])
    else:
        print " %3.3E  (%3.3E,%3.3E) & %3.3f (%3.3f,%3.3f) \\"%(KS_mean,(KS_interval)[0],(KS_interval)[1],D_mean,D_interval[0],D_interval[1])

    ax.step(np.concatenate([sorted_data_II,[0]]),np.arange(0.0,sorted_data_II.size+1)/sorted_data_II.size,\
                color='r',label="II (%s)"%len(nis_II),linewidth=2)
    ax.set_title(r"$^{56}Ni \ %s$"%ni_type.split('_')[1],size=20)
    ax.set_xlabel(r"$\mathrm{Nickel \ Mass} \ M_\odot$",size=20)
    ax.legend(loc='best')
    fig.savefig("%s_%s_hist_SESN.png"%(save_name,ni_type))
    

ax2.tick_params(labelsize=17)
fig2.subplots_adjust(left=0.08)                                                                                                                                     
fig2.subplots_adjust(right=0.95)                                                                                                                                    
fig2.subplots_adjust(top=0.97)                                                                                                                                      
fig2.subplots_adjust(bottom=0.12)                                                                                                                                      
pl.setp(pl.legend().get_texts(), fontsize='16')
ax2.legend(loc='best',prop={'size':17},ncol=1)
ax2.set_xlabel(r"$\mathrm{Nickel \ Mass} \ M_\odot$",size=22)
ax2.set_xlim(0,0.21)
fig2.savefig("%s_all_hist_SESN.png"%save_name,dpi=300)
pl.show()



for sn_type in SN_plot:
    where_type = np.where(np.logical_and(my_nis['type']==sn_type,~np.isnan(my_nis['%s'%ni_type])))[0]
    print sn_type
    for ni_type in Ni_types:
        print ni_type
        print stats.ks_2samp(my_nis["%s"%ni_type][where_type],sorted_data_II)

#    for sn,ni_peak,ni_k,ni_tail in my_nis[['SN','Mni_peak','Mni_Khatami','Mni_tail']][where_type]:
#        print sn,ni_peak,ni_k,ni_tail


where_BL = np.where(np.logical_and(my_nis['type']=="Ic_GRB",~np.isnan(my_nis['%s'%ni_type])))[0]
nis_peak_BL = my_nis[where_BL]['Mni_peak']
print np.mean([i-j for i in sorted_data_II for j in nis_peak_BL])
