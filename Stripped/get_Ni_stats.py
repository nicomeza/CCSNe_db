from astro_utils import *
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
from scipy.misc import derivative
from scipy.optimize import newton,brentq
from astroML.plotting import hist


sn_labels = ['sn','type','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_test.dat',dtype={'names':sn_labels,'formats':sn_formats})


SN_plot = {'IIb':{"marker":"o","color":"b",'string':'$\mathrm{IIb}$'},'Ib':{"marker":"p","color":"g",'string':'$\mathrm{Ib}$'},'Ic':{"marker":"s","color":"k",'string'\
:'$\mathrm{Ic}$'},'Ibc':{"marker":"x","color":"k",'string':'$\mathrm{Ibc}$'},'Ic_GRB':{"marker":"*","color":"r",'string':'$\mathrm{Ic-GRB}$'},'Ic_BL':{"marker":"s","c\
olor":"r",'string':'$\mathrm{Ic_{BL}}$'}}



colormap = pl.cm.spectral

my_file = "Ni56_BVRIYJH_peak.dat"
my_nis = np.genfromtxt(my_file,\
                           names=('SN','type','tp','Lp','logLp','Mni_peak','Mni_K','Mni_tail','trise','tdecay','FWHM','t_inf','L_dot_inf'),\
                           dtype=('S15','S10','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

nis_II = [l.split() for l in open('SNII_Ni56.dat').readlines()[1:]]
nis_II = [(sn,np.median(np.asarray(n.split(','),dtype='f8'))) for sn,n in nis_II]

where_tail = np.where(~np.isnan(my_nis['Mni_tail']))[0]
where_rise = np.where(~np.isnan(my_nis['tp']))[0]
nis_rise = my_nis[where_rise]
nis_compare = my_nis[where_tail]

fig = pl.figure()
ax = fig.add_subplot(111)

nis_line = np.arange(0.0,0.5,0.1)
ax.plot(nis_compare['Mni_peak'],nis_compare['Mni_tail'],marker='o',linestyle='None')
ax.plot(nis_line,nis_line)
ax.set_xlabel(r"$\mathrm{Arnett} \ ^{56}Ni \ [M_\odot]$",size=15)
ax.set_ylabel(r"$\mathrm{Tail} \ ^{56}Ni \ [M_\odot]$",size=15)
ax.set_xlim(0,np.max(nis_compare['Mni_peak'])+0.05)
ax.set_ylim(0,np.max(nis_compare['Mni_peak'])+0.05)
pl.show()


pl.hist(nis_compare['Mni_peak'],label='Peak',alpha=0.5)
pl.hist(nis_compare['Mni_tail'],label='Tail',alpha=0.5)
pl.legend()
pl.show()

pl.hist(nis_compare['Mni_K'],label='Khatami',alpha=0.5)
pl.hist(nis_compare['Mni_tail'],label='Tail',alpha=0.5)
pl.legend()
pl.show()

## Rise time distributions ##

for sn_type in SN_plot:
    
    where_type = np.where(nis_rise['type']==sn_type)[0]
    if len(where_type)>0:
        print sn_type, len(where_type)
        if len(where_type)>3:
            hist(nis_rise['tp'][where_type],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='knuth')
        else:
            hist(nis_rise['tp'][where_type],label="%s (%s)"%(sn_type,len(where_type)),color=SN_plot[sn_type]['color'],alpha=0.5,bins='blocks')
        pl.axvline(np.median(nis_rise['tp'][where_type]),linestyle='--',color=SN_plot[sn_type]['color'],linewidth=5)
    
pl.legend(loc='best')
pl.xlabel(r'$\mathrm{Rise \ time \ [days]}$',size=15)
pl.savefig("rise_hist_subtypes.png")
pl.show()

## Nickel distributions ####

Ni_types = ['Mni_tail','Mni_K','Mni_peak'] 

for ni_type in Ni_types:

    fig = pl.figure()
    ax = fig.add_subplot(111)

    for sn_type in SN_plot:
        
        #if sn_type not in ['Ic_BL','Ic_GRB','Ibc','Ib']:
        where_type = np.where(np.logical_and(my_nis['type']==sn_type,~np.isnan(my_nis['%s'%ni_type])))[0]
        if len(where_type)>1:
            sorted_data = np.sort(my_nis['%s'%ni_type][where_type]) 
            ax.step(np.concatenate([sorted_data,[0]]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,color=SN_plot[sn_type]['color'],label="%s (%s)"%(sn_type,len(where_type)))
            #ax.hist(my_nis['%s'%ni_type][where_type],color=SN_plot[sn_type]['color'],label=sn_type,alpha=0.1,ls='dashed',cumulative=True,normed=True,fill=None)

    pl.title(ni_type)
    pl.legend(loc='upper left')
    pl.savefig("%s_hist_subtypes.png"%(ni_type))
    pl.show()

    
    
for ni_type in Ni_types:

    fig = pl.figure()
    ax = fig.add_subplot(111)
   
    sn_type = "IIb"
    where_type = np.where(np.logical_and(my_nis['type']==sn_type,~np.isnan(my_nis['%s'%ni_type])))[0]
    if len(where_type)>1:
        print "Number of %s : %s"%(sn_type,len(where_type))
        sorted_data = np.sort(my_nis['%s'%ni_type][where_type]) 
        ax.step(np.concatenate([sorted_data,[0]]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,\
                    color=SN_plot[sn_type]['color'],label="IIb (%s)"%(len(where_type)))
    
    sn_type = "Ic"
    where_type = np.where(np.logical_and(np.logical_and(my_nis['type']!="IIb",my_nis['type']!="Ic_GRB"),~np.isnan(my_nis['%s'%ni_type])))[0]
    if len(where_type)>1:
        print "Number of %s : %s"%(sn_type,len(where_type))
        sorted_data = np.sort(my_nis['%s'%ni_type][where_type]) 
        ax.step(np.concatenate([sorted_data,[0]]),np.arange(0.0,sorted_data.size+1)/sorted_data.size,\
                    color=SN_plot[sn_type]['color'],label="Ibc (%s)"%len(where_type))
        

    sorted_data_II = np.sort([x[1] for x in nis_II]) 
    ax.step(np.concatenate([sorted_data_II,[0]]),np.arange(0.0,sorted_data_II.size+1)/sorted_data_II.size,\
                color='r',label="II (%s)"%len(nis_II))
    
    pl.title(ni_type)
    pl.legend(loc='upper left')
    ax.set_xlabel(r"$\mathrm{Nickel \ Mass}$")
    pl.savefig("%s_hist_IIb-Ibc.png"%(ni_type))
    pl.show()



for sn_type in SN_plot:
    where_type = np.where(np.logical_and(my_nis['type']==sn_type,~np.isnan(my_nis['%s'%ni_type])))[0]
    print sn_type
    for sn,ni_peak,ni_k,ni_tail in my_nis[['SN','Mni_peak','Mni_K','Mni_tail']][where_type]:
        print sn,ni_peak,ni_k,ni_tail
