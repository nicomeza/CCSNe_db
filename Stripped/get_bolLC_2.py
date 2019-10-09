# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pl
from astro_utils import *
import os
import argparse
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
from scipy.stats import entropy
import bisect 

def get_lc(filter_file,use_filters,t_0,t_non,t_discov,t_upper=90.):
    
    scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = None,None,None

    if os.path.isfile(filter_file):
        print filter_file,filter
        try:
            scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = np.genfromtxt(filter_file,usecols=[0,1,2]).T
            if t_0<1000 and t_0>0.0:
                if t_non>1000 and t_discov>1000:
                    t_0 = 0.5*(t_non + t_discov)
                    print "Taking t_0 as midpoint between last non detection and discovery"
                    no_nebular = np.where(np.logical_and(scope['jd_%s'%filter]<(t_0+t_upper),scope['err_%s'%filter]<0.35))[0]
                    if len(no_nebular)>1:
                        
                        scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = \
                            scope['jd_%s'%filter][no_nebular]-t_0,scope['mag_%s'%filter][no_nebular],scope['err_%s'%filter][no_nebular] 
                        sort_arg = np.argsort(scope['jd_%s'%filter])
                        scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = \
                            scope['jd_%s'%filter][sort_arg],scope['mag_%s'%filter][sort_arg],scope['err_%s'%filter][sort_arg]
                    else:
                        print scope['jd_%s'%filter][no_nebular]-t_0
                        filter_removed = use_filters.pop(use_filters.index(filter))
                        print "LACK OF DATA : filter %s removed"%(filter_removed)

                else:
                    print "Cant use. No t_0 and no discovery/non detection"
                    filter_removed = use_filters.pop(use_filters.index(filter))
                    print "LACK OF DATA : filter %s removed"%(filter_removed)
                    
            else:    
                no_nebular = np.where(np.logical_and(scope['jd_%s'%filter]<(t_0+t_upper),scope['err_%s'%filter]<0.35))[0]
                
                if len(no_nebular)>1:
                    
                    scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = \
                        scope['jd_%s'%filter][no_nebular]-t_0,scope['mag_%s'%filter][no_nebular],scope['err_%s'%filter][no_nebular] 
                    sort_arg = np.argsort(scope['jd_%s'%filter])
                    scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = \
                        scope['jd_%s'%filter][sort_arg],scope['mag_%s'%filter][sort_arg],scope['err_%s'%filter][sort_arg]
                    
                else:
                    print scope['jd_%s'%filter][no_nebular]-t_0
                    filter_removed = use_filters.pop(use_filters.index(filter))
                    print "LACK OF DATA : filter %s removed"%(filter_removed)
                    
                    
        except:
            
            print "No errors in band %s ?"%filter
            filter_removed = use_filters.pop(use_filters.index(filter))
            print "filter %s removed"%(filter_removed)
            scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = np.genfromtxt(filter_file,usecols=[0,1,2]).T
            no_nebular = np.where(scope['jd_%s'%filter]<(t_0+t_upper))[0]
            if len(no_nebular)>2:
                scope['jd_%s'%filter],scope['mag_%s'%filter] = scope['jd_%s'%filter][no_nebular]-t_0,scope['mag_%s'%filter][no_nebular]
            
    else:
        print "no %s"%filter_file
        filter_removed = use_filters.pop(use_filters.index(filter))
        print "filter %s removed"%(filter_removed)

    return scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter]


def get_boundaries(FILTER_ZPS,show=False,plot=True,SN_name="None"):
    
    tmin_max,tmax_min = 10000. , -1.0
    max_cadence,min_cadence = 0.,100.0
    min_cadence_filt,max_cadence_filt = None,None
    
    if plot:
        pl.gca().set_color_cycle(filters_colors)

    filters_cadence_full = []
    offset = 0.0
    
    for filter in FILTER_ZPS:
        
        lco_filter = filter['filter']
        lco_name = lco_filter#lcogt_names[lcogt_filters.index(lco_filter)]
        print lco_filter
        print filter['lambda']
                
        try:
            scope['mag_%s'%lco_name]
            if len(scope['mag_%s'%lco_name])<=1:
                continue
        except:
            continue
        
        mag = scope['mag_%s'%lco_name] 
        err = scope['err_%s'%lco_name]
        t = scope['jd_%s'%lco_name]

        if plot:
            split_filter = lco_filter.split('_')
            if len(split_filter)>1:
                if split_filter[1] in ['CSPI','3014','3009','9844']:
                    pl.errorbar(t,mag+offset,yerr=err,fmt='o',linestyle='None',label="%s (%s)"%(split_filter[0],"Swope"))
                elif split_filter[1] in ['WIRC']:
                    pl.errorbar(t,mag+offset,yerr=err,fmt='o',linestyle='None',label="%s (%s)"%(split_filter[0],"du Pont"))
                elif split_filter[1] in ['RC','RC1','RC2']:
                    pl.errorbar(t,mag+offset,yerr=err,fmt='o',linestyle='None',label="%s (%s)"%(split_filter[0],"RetroCam"))
                else:
                    pl.errorbar(t,mag+offset,yerr=err,fmt='o',linestyle='None',label="%s (%s)"%(split_filter[0],split_filter[1]))
            else:
                pl.errorbar(t,mag+offset,yerr=err,fmt='o',linestyle='None',label=lco_filter)
            offset += 1.0

        tmax,tmin = np.max(t),np.min(t)
        dtmax = tmax-tmin
        if tmin >= tmax_min:
            tmax_min = tmin
        if tmax <= tmin_max:
            tmin_max = tmax
            
        scope['min_%s'%lco_filter] ,scope['max_%s'%lco_filter]  = tmin,tmax
        scope['cadence_%s'%lco_filter] = len(t)/dtmax   #NEED CADENCE IN THE ACTUAL INTERVAL (tmax_min,tmin_max)
            #scope['inter_%s'%lco_filter] = interp1d(t,mag,kind='linear')
            #pl.plot(t,scope['inter_%s'%lco_filter](t),color='k',linestyle='--')
        
        filters_cadence_full.append(len(t)/dtmax)
        print len(t)/dtmax, (tmin,tmax)
        if max_cadence < len(t)/dtmax:
            max_cadence_filt = lco_name
            max_cadence = len(t)/dtmax
            
        if min_cadence > len(t)/dtmax:
            min_cadence_filt = lco_name
            min_cadence = len(t)/dtmax
    
    if plot:

        pl.axvline(tmin_max,linestyle='--',color='k')
        pl.axvline(tmax_min,linestyle='--',color='k')
        pl.xlabel(r'$\mathrm{Time \ since \ explosion}$')
        pl.ylabel(r'$\mathrm{Magnitude \ (arbitrary)}$')
        pl.title(r'$%s$'%(SN_name))
        pl.legend(loc='best',ncol=2,prop={'size':9})
        pl.gca().invert_yaxis()
        pl.savefig("%s_lcs.png"%SN)
    
        if show:
            pl.show()
        else:
            pl.close()

    max_cadence,min_cadence = 0.,100.0
    min_cadence_filt,max_cadence_filt = None,None
    filters_cadence = []
    filters_entropy = []

    for filter in FILTER_ZPS:
            
        lco_filter = filter['filter']
        lco_name = lco_filter
        try:
            scope['mag_%s'%lco_name]
            if len(scope['mag_%s'%lco_name])<=1:
                continue
        except:
            continue
        
        t = scope['jd_%s'%lco_name]
        where_cadence = np.where(np.logical_and(t <= tmin_max + 0.8 ,t >= tmax_min - 0.8))[0]
        if len(where_cadence)<2:
            print "No data in range (tmax_min,tmin_max) in filter %s"%filter
            continue
        else:
            dtmax = np.max(t[where_cadence]) - np.min(t[where_cadence])
            scope['cadence_%s'%lco_filter] = len(where_cadence)/dtmax   #NEED CADENCE IN THE ACTUAL INTERVAL (tmax_min,tmin_max)
            filter_cadence = len(where_cadence)/dtmax
            filters_cadence.append(filter_cadence)
            n,bins = np.histogram(t[where_cadence])
            filters_entropy.append(entropy(n))
            print filter_cadence
            if max_cadence < filter_cadence:
                max_cadence_filt = lco_name
                max_cadence = filter_cadence
                
            if min_cadence > filter_cadence:
                min_cadence_filt = lco_name
                min_cadence = filter_cadence
        
    print "Max cadence filter between (%s,%s) : %s  "%(tmax_min,tmin_max,max_cadence_filt)
    print "Min cadence filter between (%s,%s) : %s  "%(tmax_min,tmin_max,min_cadence_filt)
    max_entropy_filter =  FILTER_ZPS[np.argmax(filters_entropy)]['filter']
    print "Max entropy filter between (%s,%s) : %s  "%(tmax_min,tmin_max,max_entropy_filter)

    return tmax_min,tmin_max,max_entropy_filter


def interp_lc(where_to_t,t,data,dataerr,inter_t,inter_data,inter_dataerr):

    try:
        for t_inter in where_to_t:
                
            nearest_t = np.where(np.abs(t-t_inter)<0.25)[0]
        
            if len(nearest_t)>0:
                nearest_t = np.argmin(np.abs(t-t_inter)) 
                inter_t.append(t[nearest_t])
                inter_data.append(data[nearest_t])
                inter_dataerr.append(dataerr[nearest_t])
                
            else:        
                inter_t.append(t_inter)
                inter,inter_err = inter_sample(t_inter,t,data,yerr=dataerr,n_sample=500,plot=False,confidence=False,honest_mean=True)
                inter_data.append(inter)
                inter_dataerr.append(inter_err)
        return True
    except:
        print "Something went wrong in interpolation."
        return False


L_sun = 3.8270e33
M_sun = 4.7554
Ni_decay = 8.8 # days                                                                                                   
CO_decay = 111.3 # days                                                
R_sun = 6.96e10 # cm                                                                
pc_to_cm = 3.086e+18
day_to_sec = 60*60*24.0
M_sun_gr = 1.989e33 # Sun mass in grams                                                              
c_kms = 299792.458 #km/s                                                                                  
c_cms = 2.99792458*1e10
pi = np.pi
R_V = 3.1
R_U = 1.569*R_V # From Cardelli law                                                                                            
R_B = 4.2
R_R = 0.75*R_V
R_I = 0.478*R_V
R_u = R_V*1.579
R_g = R_V*1.161
R_r = R_V*0.843
R_i = R_V*0.639
R_z = R_V*0.453  # from shlegel et al (see Kim-Lee 2007)                            

scope = globals()
parser = argparse.ArgumentParser(description='Calculate bolometric light curves.')
parser.add_argument('Args', metavar='Img', type=str, nargs='+',help='--')
args = parser.parse_args()
IN_FILE = args.Args[0] # SN metadata

# TO USE IN IPYTHON : 
# import sys
# sys.argv = ['arg1','arg2']


#------------------ AB ZERO POINT -------------------------------------

AB_ZPs = np.genfromtxt('./AB_zeropoints.dat',usecols=[0,1,2,3,4],dtype={'names':('filter','lambda','lambda_p','flux','ZP'),'formats':('S10','f8','f8','f8','f8')})

AB_ZPs['lambda']*=1e4 # Transforming wavelength to Angstroms


#---------- CSP zeropoints ------------------------------

CSP_ZPs = np.genfromtxt('./CSP/CSP_AB_ZPs.dat',dtype={'names':('filter','ZP'),'formats':('S10','f8')})

#---------- CSP SN data-----------------------------------

csp_names = ('SN', 'ra', 'dec', 'host', 'host_type', 'host_z_hel', 'discovery', 'discoverer', 'SN_type', 'phase_i', 'phase_f' ,'maximum_epoch', 'maximum_epoch_err')
csp_types = ('S10','S20','S20','S30','S20','f8','S20','S20','S6','S5','S5','f8','f8')
CSP_SN = np.genfromtxt('./CSP/CSP_SN.dat',usecols=np.arange(len(csp_names)),dtype={'names':csp_names,'formats':csp_types})

# -------------SN METADATA --------------------

show = False
host_redd = True

print " Convert photometry to fluxes to obtain SED "


sn_labels = ['sn','type','host','hostredshift','hostlumdist','mw_ebv','host_ebv','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8']

    
SN_DATA = np.genfromtxt(IN_FILE,dtype={'names':sn_labels,'formats':sn_formats})

current_dir = os.getcwd()
print "Current directory : %s"%current_dir
colormap = pl.cm.spectral

filters = ['W2_uvot','M2_uvot','W1_uvot','U_uvot','B_uvot','V_uvot','U','B','V','R','I','u','g','r','i','z','J','H','K','Ks']
filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(filters))]
        
#use_filters = ['W2_uvot','M2_uvot','W1_uvot','U','B','V','R','I','J','H','K','Ks']
#use_filters = ['W2_uvot','M2_uvot','W1_uvot','U_uvot','B','V','R','I','J','H','K','Ks']

Ni_56 = []
Lps = []
tps = []


for SN,z_SN,mw_E_B_V,host_E_B_V,t_0,d_L,t_discov,t_non in SN_DATA[['sn','sn_z','mw_ebv','host_ebv','t_0','hostlumdist','t_discov','t_non_det']]:

    print "####### %s ########## \n"%SN
    use_filters = ['B_CSPI','B','B_AB','V','V_3014','V_3009','V_9844','V_AB','R','R_AB','I','I_AB','g','g_CSPI','g_AB','r','r_CSPI','r_AB','i','i_CSPI','i_AB','R_c',\
                       'J','H','J_AB','H_AB','Y_WIRC','Y_RC','J_WIRC','H_WIRC','J_RC1','J_RC2','Y_RC','H_RC']
    #use_filters = ['B','V','R','I','R_c']
    #use_filters = ['B_uvot','V_uvot','R','I','R_c']
    try: 
        sn_log = open("%s.log"%SN,'w')
        os.chdir("%s/"%SN)
        
        # ---------------------LOAD PHOTOMETRY ----------------------------------------------------

        for filter in zip(use_filters):
            filter = filter[0]
            filter_file = "mags_%s.out"%filter
            scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = get_lc(filter_file,use_filters,t_0,t_non,t_discov)
            
        FILTER_ZPS = AB_ZPs[np.in1d(AB_ZPs['filter'],use_filters)]  
        print FILTER_ZPS['filter'],use_filters

        if len(use_filters)<=1:
            print "NOT ENOUGH FILTERS AVAILABLE"
            os.chdir(current_dir)
            continue
        filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(FILTER_ZPS))]
        
        print "# ---------- FILTER_ZPS filters ----------------------------"

        tmax_min,tmin_max,max_entropy_filter = get_boundaries(FILTER_ZPS,SN_name=SN,plot=False)
        
        baseline = np.arange(int(tmax_min)+1,int(tmin_max),1)
        #inter_t = scope['jd_%s'%min_cadence_filt]
        inter_t = scope['jd_%s'%max_entropy_filter]
        max_t_inter,min_t_inter = np.max(inter_t),np.min(inter_t)
        where_inter = np.where(np.logical_and(inter_t <= tmin_max + 0.8 ,inter_t >= tmax_min - 0.8))[0]

        if max_t_inter>tmin_max: 
            next_inter = inter_t[where_inter[-1]+1]
            if (next_inter - inter_t[where_inter[-1]])<15.:
                new_mag,new_err = inter_sample(tmin_max,inter_t,scope['mag_%s'%max_entropy_filter],yerr=scope['err_%s'%max_entropy_filter],\
                                                            n_sample=500,plot=False,confidence=False,honest_mean=True)
                print "Inserting (%2.2f,%2.2f,%2.2f) in %s"%(tmin_max,new_mag,new_err,max_entropy_filter)
                scope['jd_%s'%max_entropy_filter] = list(scope['jd_%s'%max_entropy_filter])
                scope['mag_%s'%max_entropy_filter] = list(scope['mag_%s'%max_entropy_filter])
                scope['err_%s'%max_entropy_filter] = list(scope['err_%s'%max_entropy_filter])
                bisect.insort(scope['jd_%s'%max_entropy_filter],tmin_max) 
                insert_index = scope['jd_%s'%max_entropy_filter].index(tmin_max)
                scope['jd_%s'%max_entropy_filter] = np.array(scope['jd_%s'%max_entropy_filter])
                scope['mag_%s'%max_entropy_filter] = np.array(scope['mag_%s'%max_entropy_filter][0:insert_index]+[new_mag]+scope['mag_%s'%max_entropy_filter][insert_index:])
                scope['err_%s'%max_entropy_filter] = np.array(scope['err_%s'%max_entropy_filter][0:insert_index]+[new_err]+scope['err_%s'%max_entropy_filter][insert_index:])

        inter_t = scope['jd_%s'%max_entropy_filter]
        where_inter = np.where(np.logical_and(inter_t <= tmin_max + 0.8 ,inter_t >= tmax_min - 0.8))[0]
        
        if min_t_inter<tmax_min: 
            
            prev_inter = inter_t[where_inter[0]-1]
            if (-prev_inter + inter_t[where_inter[0]])<10.:
                new_mag,new_err = inter_sample(tmax_min,inter_t,scope['mag_%s'%max_entropy_filter],yerr=scope['err_%s'%max_entropy_filter],\
                                                            n_sample=500,plot=False,confidence=False,honest_mean=True)
                print "Inserting (%2.2f,%2.2f,%2.2f) in %s"%(tmax_min,new_mag,new_err,max_entropy_filter)
                scope['jd_%s'%max_entropy_filter] = list(scope['jd_%s'%max_entropy_filter])
                scope['mag_%s'%max_entropy_filter] = list(scope['mag_%s'%max_entropy_filter])
                scope['err_%s'%max_entropy_filter] = list(scope['err_%s'%max_entropy_filter])
                bisect.insort(scope['jd_%s'%max_entropy_filter],tmax_min) 
                insert_index = scope['jd_%s'%max_entropy_filter].index(tmax_min)
                scope['jd_%s'%max_entropy_filter] = np.array(scope['jd_%s'%max_entropy_filter])
                scope['mag_%s'%max_entropy_filter] = np.array(scope['mag_%s'%max_entropy_filter][0:insert_index]+[new_mag]+scope['mag_%s'%max_entropy_filter][insert_index:])
                scope['err_%s'%max_entropy_filter] = np.array(scope['err_%s'%max_entropy_filter][0:insert_index]+[new_err]+scope['err_%s'%max_entropy_filter][insert_index:])
                
        inter_t = scope['jd_%s'%max_entropy_filter]
        where_inter = np.where(np.logical_and(inter_t <= tmin_max + 0.8 ,inter_t >= tmax_min - 0.8))[0]
        
        sn_log.write("Filters used:")
        sn_log.write("\t".join(FILTER_ZPS['filter'])+"\n")
        sn_log.write("Max entropy filter : %s \n"%max_entropy_filter)
        sn_log.write("Number of data points : %s \n"%len(where_inter))
        sn_log.write("Data points time range : (%3.3f,%3.3f) \n"%(tmax_min,tmin_max))
        if len(where_inter)<2:
            print "No data in range (tmax_min,tmin_max)"
            sn_log.write("No data in range (tmax_min,tmin_max) \n")
            os.chdir(current_dir)
            sn_log.close()
            continue
        if tmax_min>=tmin_max:
            print "tmax_min>=tmin_max. Bad overlap of filters."
            sn_log.write("tmax_min>=tmin_max. Bad overlap of filters.\n")
            os.chdir(current_dir)
            sn_log.close()
            continue

        
        print "# ---------- filter interpolation : ----------------------------"

        filters_used = []
        
        for filter in FILTER_ZPS : 
    
            lco_filter = filter['filter']
            lco_name = lco_filter #lcogt_names[lcogt_filters.index(lco_filter)]
            x0 = 1/(filter['lambda']*1e-4)
            R_x = (R_V*a_x(x0)+b_x(x0))
    
            try:
                if host_E_B_V<99 and host_redd:
                    mag = scope['mag_%s'%lco_name] - R_x*(mw_E_B_V+host_E_B_V)   # Deredenned magnitudes 
                else:
                    mag = scope['mag_%s'%lco_name] - R_x*mw_E_B_V   # Deredenned magnitudes 
            except:
                continue
                print "Fail"
                                
            err = scope['err_%s'%lco_name]
            t = scope['jd_%s'%lco_name]
            scope['inter_t_%s'%lco_name] = []
            scope['inter_mags_%s'%lco_name] = []
            scope['inter_err_%s'%lco_name] = []
            filters_used.append(lco_name)
            
            interp_boolean = interp_lc(inter_t[where_inter],t,mag,err,scope['inter_t_%s'%lco_name],scope['inter_mags_%s'%lco_name],scope['inter_err_%s'%lco_name])

            if interp_boolean:
                print "%s interpolated"%lco_filter
            else:
                print "%s interpolation failed"%lco_filter
            
            
        inter_mags = [ (filter,np.array(scope['inter_mags_%s'%filter]),np.array(scope['inter_err_%s'%filter])) \
                   for filter in filters_used]

        sorted_mags = sorted(inter_mags,\
                                key = lambda x : FILTER_ZPS['lambda'][np.where(FILTER_ZPS['filter']==x[0])[0]])

        sorted_lambdas = np.sort(FILTER_ZPS['lambda_p'][np.in1d(FILTER_ZPS['filter'],filters_used)])
        sorted_bands = [FILTER_ZPS['filter'][np.where(FILTER_ZPS['lambda_p']==lam)[0]][0] for lam in sorted_lambdas]

        fluxes = [[] for i in range(len(where_inter))]
        fluxes_err = [[] for i in range(len(where_inter))]
        
        filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(sorted_mags))]
        pl.gca().set_color_cycle(filters_colors)
        
        print "# ---------- Getting fluxes ----------------------------"
        
        for filter,mags,errs in sorted_mags:
            
            print filter
            index = filters_used.index(filter)
            lambda_f = FILTER_ZPS['lambda_p'][index]
            ZP = FILTER_ZPS['ZP'][index]
            print lambda_f,ZP
            
            # m_AB = -2.5log(f_AB) - 48.6 - ZP
            F_AB = 10**(-0.4*(48.6+mags+ZP)) # Jansky 
            F_err =  F_AB*0.4*errs*np.log(10)
            
            pl.errorbar(inter_t[where_inter],mags,yerr = errs,label=filter,linestyle='None',fmt='o')
    
            for i,f in enumerate(zip(F_AB,F_err)):
                
                fab = f[0]
                ferr = f[1]
                fluxes[i].append(fab)
                fluxes_err[i].append(np.asarray(ferr))
                
        pl.axvline(tmin_max,linestyle='--',color='k')
        pl.axvline(tmax_min,linestyle='--',color='k')
        pl.legend(loc='best',ncol=2,prop={'size':8})
        pl.gca().invert_yaxis()
        pl.savefig("%s_inter_lcs.png"%SN)
        if show:
            pl.show()
        else:
            pl.close()
        
        
        print "# --------- Plot SED's --------------------------------------"
    

        fig = pl.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        
        colormap = pl.cm.spectral
        fig.gca().set_color_cycle([colormap(i) for i in np.linspace(0.1, 0.9,len(where_inter))])
        
        filter_names = AB_ZPs['filter'][np.in1d(AB_ZPs['lambda_p'],sorted_lambdas)]
        
        for epoch,flux,flux_err in zip(inter_t[where_inter],fluxes,fluxes_err):
            
            flux = np.asarray(flux)
            ax.errorbar(c_cms/(sorted_lambdas*1e-8),flux,yerr=flux_err,linestyle='--',label='%2.2f'%epoch)
                     
        ax.set_title(SN)
        ax.set_xlabel('Frequency [Hz]',size=22)
        ax.set_ylabel('Flux in ergs/s/cm2/Hz',size=22)
        ax.legend(loc='best',ncol=2,prop={'size':10})
        fig.savefig('./SED_nu_%s.png'%R_V)

        if show:
            pl.show()
        else:
            pl.close()

        fig = pl.figure(figsize=(10,8))
        ax = fig.add_subplot(111)

        colormap = pl.cm.spectral
        fig.gca().set_color_cycle([colormap(i) for i in np.linspace(0.1, 0.9, len(where_inter))])
        prev_epoch = -10
        for epoch,flux,flux_err in zip(inter_t[where_inter],fluxes,fluxes_err):
            if (epoch-prev_epoch)>0.5:  
                flux = np.asarray(flux)
                ax.errorbar(sorted_lambdas,flux*c_cms/((sorted_lambdas*1e-4)**2.),\
                                yerr=c_cms*np.asarray(flux_err)/((sorted_lambdas*1e-4)**2.),marker='o',linestyle='--',label='%2.1f'%epoch,alpha=0.7)
                prev_epoch = epoch

        trans = ax.get_xaxis_transform()

        band_string = ""
        for lam,flux in zip(sorted_lambdas,fluxes[0]):
            flux = np.asarray(flux)
            max_flux = np.max(flux)*c_cms/((lam*1e-4)**2.)
            band = FILTER_ZPS['filter'][np.where(FILTER_ZPS['lambda_p']==lam)[0]]
            band_string+=str(band[0].split('_')[0])
            print band,lam,max_flux
            ax.axvline(lam,linestyle='--',color='k',alpha=0.4)
            ax.annotate(r"$\mathrm{%s}$"%str(band[0].split("_")[0]),(lam+5,1.03),xycoords=trans,size=21,rotation=0)
            
        max_flux = np.max(np.array(fluxes)*c_cms/((sorted_lambdas*1e-4)**2.))
        print "upper limit = %s"%max_flux
        #ax.set_title(r"$\mathrm{%s}$"%SN,size=17)
        ax.annotate(r"$\mathrm{%s}$"%SN,(8100,0.95),xycoords=trans,size=22)
        ax.set_xlabel(r'$\mathrm{Rest \ Wavelength \ [\AA]}$',size=22)
        ax.set_ylabel(r'$\mathrm{Flux \ in \ ergs/s/cm^2/\AA}$',size=22)
        ax.set_xlim(4000,16500)
        ax.set_ylim(0,1.05*max_flux)
        ax.legend(loc='best',ncol=2,prop={'size':15})
        fig.subplots_adjust(left=0.07,right=0.97,bottom=0.1,top=0.94)
        ax.tick_params(labelsize=17)
        fig.savefig('./SED_lam_%s_%s.png'%(band_string,str(R_V)),dpi=300)
        if show:
            pl.show()
        pl.close()
    
        print "# --------- Saving SEDs -------------------------------------"
        
        d_SED = {}
        d_SED['t'] = inter_t[where_inter]
        d_SED['data'] = fluxes
        d_SED['daterr'] = fluxes_err
        d_SED['lambda'] = sorted_lambdas
        d_SED['bands'] = sorted_bands 
        np.savez( 'SED_%s_%s.npz'%(SN,band_string),**d_SED)
        
        
        print "------------ Integrating SED -------------------------"
        
        fig = pl.figure(figsize=(10,6))
        ax = fig.add_subplot(111)

        from scipy.integrate import simps
        L_bol = []
        M_bol = []
        Lfile = open('%s_Lbol_%s.dat'%(SN,band_string),'w')
        if host_redd:
            Lfile.write("## SN \t z \t MW_E_B_V \t host E_B_V d [Mpc] \t t_0 \n")
            Lfile.write("## %s \t %s \t %s \t  %s \t %s \t %s\n"%(SN,z_SN,mw_E_B_V,host_E_B_V,d_L,t_0))
        else:
            Lfile.write("## SN \t z \t MW_E_B_V \t d [Mpc] \t t_0 \n")
            Lfile.write("## %s \t %s \t  %s \t %s \t %s\n"%(SN,z_SN,mw_E_B_V,d_L,t_0))
        Lfile.write("###############%s#######################\n"%band_string)
        Lfile.write("# t-t_0 \t Lbol [erg/s] \t log(Lbol) \t M_ni\n")
        int_nu = False
        
        M_ni = lambda t_r,L_bol : L_bol/Q_t(t_r,1)
        
        for epoch,flux,flux_err in zip(inter_t[where_inter],fluxes,fluxes_err):

            if int_nu:
                y,x = flux,c_cms/(sorted_lambdas*1e-8)
                I = -np.trapz(y,x)
                L = 4*pi * I * (d_L*1e6*pc_to_cm)**2.
                print epoch,I,4*pi*I*(d_L*1e6*pc_to_cm)**2. 
                L_bol.append(L)
                ax.plot(epoch,L,marker='o',color='k')
                Lfile.write("%s \t %s \n"%(epoch,L))
                
            else:

                y,x = (1+z_SN)*np.asarray(flux)*c_cms/(sorted_lambdas*1e-4)**2.,sorted_lambdas
                I = simps(y,x)
                L = 4*pi * I * (d_L*1e6*pc_to_cm)**2.
                print epoch,I,L,M_ni(epoch,L) 
                L_bol.append(L)
                ax.plot(epoch,L,marker='o',color='k')
                Lfile.write("%s \t %s \t %s \t %s\n"%(epoch,L,np.log10(L),M_ni(epoch,L)))
                
        L_bol = np.asarray(L_bol)
        where_2nd_peak = np.where(inter_t[where_inter]>10.)[0]
        where_peak = np.argmax(L_bol[where_2nd_peak])
        try:
            in_peak = np.where(L_bol==L_bol[where_2nd_peak][where_peak])[0][0]
            sn_log.write("Peak at %2.2f days past explosion.\n"%inter_t[where_inter][where_2nd_peak][where_peak])
            print "Number of data points before/after apparent peak (after 10 days past t_0) : (%s/%s)"%(len(L_bol[:in_peak]),len(L_bol[in_peak:]))
            sn_log.write("Number of data points before/after apparent peak (after 10 days past t_0) : (%s/%s)"%(len(L_bol[:in_peak]),len(L_bol[in_peak:])))
        except:
            print "Could not get number of data points before/after peak."
        tp,Lp = inter_t[where_inter][where_2nd_peak][where_peak],L_bol[where_2nd_peak][where_peak]
        M_p = M_ni(tp,Lp)
        print "t_peak = %s\nL_peak = %s\nM_ni = %s\n"%(tp,Lp,M_p)
        
        Lfile.write("#t_peak = %s\n#L_peak = %s\n#M_ni = %s\n"%(tp,Lp,M_p))
        Lfile.close()
        ax.set_xlabel('t')
        ax.set_ylabel('L  [erg/s]')
        fig.savefig('./Lbol_%s_%s.png'%(SN,band_string))
        if show:
            pl.show()
        pl.close()
        
        Ni_56.append(M_p)
        Lps.append(Lp)
        tps.append(tp)
        
        try:
            os.chdir(current_dir)
            sn_log.close()
            pl.close()
        except:
            sn_log.close()
            pl.close()
            print "No diretory %s"%current_dir

    except:
        
        print "No directory %s/"%SN
        os.chdir(current_dir)
        sn_log.close()
        pl.close()
        continue

try: 

    os.chdir(current_dir)
    
except:

    print "No diretory %s"%current_dir


