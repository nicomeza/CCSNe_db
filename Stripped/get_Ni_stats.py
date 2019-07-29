from astro_utils import *
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
from scipy.misc import derivative
from scipy.optimize import newton,brentq

sn_labels = ['sn','type','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_test.dat',dtype={'names':sn_labels,'formats':sn_formats})


SN_plot = {'IIb':{"marker":"o","color":"b",'string':'$\mathrm{IIb}$'},'Ib':{"marker":"p","color":"g",'string':'$\mathrm{Ib}$'},'Ic':{"marker":"s","color":"k",'string'\
:'$\mathrm{Ic}$'},'Ibc':{"marker":"x","color":"k",'string':'$\mathrm{Ibc}$'},'Ic_GRB':{"marker":"*","color":"r",'string':'$\mathrm{Ic-GRB}$'},'Ic_BL':{"marker":"s","c\
olor":"r",'string':'$\mathrm{Ic_{BL}}$'}}



colormap = pl.cm.spectral

my_file = "Ni56_BVRIYJH_peak.dat"
my_nis = np.genfromtxt(my_file,\
                           names=('SN','tp','Lp','logLp','Mni','Mni_K','Mni_tail','trise','tdecay','FWHM','t_inf','L_dot_inf'),\
                           dtype=('S15','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))


where_tail = np.where(~np.isnan(my_nis['Mni_tail']))[0]
nis_compare = my_nis[where_tail]

fig = pl.figure()
ax = fig.add_subplot(111)

nis_line = np.arange(0.0,0.5,0.1)
ax.plot(nis_compare['Mni'],nis_compare['Mni_tail'],marker='o',linestyle='None')
ax.plot(nis_line,nis_line)
ax.set_xlabel(r"$\mathrm{Arnett} \ ^{56}Ni \ [M_\odot]$",size=15)
ax.set_ylabel(r"$\mathrm{Tail} \ ^{56}Ni \ [M_\odot]$",size=15)
ax.set_xlim(0,np.max(nis_compare['Mni'])+0.05)
ax.set_ylim(0,np.max(nis_compare['Mni'])+0.05)
pl.show()


pl.hist(nis_compare['Mni'],label='Peak',alpha=0.7)
pl.hist(nis_compare['Mni_tail'],label='Tail',alpha=0.7)
pl.legend()
pl.show()

pl.hist(nis_compare['Mni_K'],label='Khatami',alpha=0.7)
pl.hist(nis_compare['Mni_tail'],label='Tail',alpha=0.7)
pl.legend()
pl.show()

