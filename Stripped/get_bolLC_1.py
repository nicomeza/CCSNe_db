# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pl
from astro_utils import *
import os

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

#--------------DECAM FILTERS-----------------------------------------------  

f_name=('lambda','u','g','r','i','z','Y','atm')
DECAM_filters = ('u','g','r','i','z','Y')                                                             
f_format = ('f8','f8','f8','f8','f8','f8','f8','S8')
DECAM = np.genfromtxt('/home/dustspeck/SN/DECAM_Filters.dat',dtype={'names':f_name,'formats':f_format})
DECAM['lambda']*=10

# ------------- LCOGT FILTERS ----------------------------------------------

LCOGT_filters = ('B','V','gp','rp','ip')
LCOGT_filters = ('gp','rp','ip')

LCOGT_DIR = "/home/dustspeck/SN/lcogt_filters"
f_name = ('lambda','response')
f_format = ('f8','f8')
lcogt_responses = []

for filter_name in LCOGT_filters:
    
    scope['%s_response'%filter_name] = np.genfromtxt("%s/SDSS.%s.txt"%(LCOGT_DIR,filter_name),dtype={'names':f_name,'formats':f_format})
    s_x,s_y = (10. * scope['%s_response'%filter_name]['lambda'],scope['%s_response'%filter_name]['response'])
    where = np.where(s_y>0)[0]
    lcogt_responses.append((s_x[where]*1e-8,s_y[where]))
    
f_name = ('filter','lambda','fwhm')
f_format = ('S3','f8','f8')   

LCOGT_data=  np.genfromtxt("%s/LCOGT_filters.dat"%(LCOGT_DIR),dtype={'names':f_name,'formats':f_format})

# ----------- UVOT FILTERS -----------------------------------------------

UVOT = ('UVW2','UVM2','UVW1','U','B','V')
UVOT_names = ('W2','M2','W1','U','B','V')
UVOT_DIR = "/home/dustspeck/SN/SWIFT/filters"
f_name = ('lambda','response')
f_format = ('f8','f8')
uvot_responses = []

for uvot_filter,uvot_name in zip(UVOT,UVOT_names):
    
    scope['%s_response'%uvot_name] = np.genfromtxt("%s/Swift_UVOT.%s.dat"%(UVOT_DIR,uvot_filter),dtype={'names':f_name,'formats':f_format})
    s_x,s_y = (scope['%s_response'%uvot_name]['lambda'],scope['%s_response'%uvot_name]['response'])
    where = np.where(s_y>0)[0]
    uvot_responses.append((s_x[where]*1e-8,s_y[where]))

#------------------ UVOT ZERO POINT -------------------------------------

AB_ZPs = np.genfromtxt('/home/dustspeck/SN/AB_zeropoints.dat',usecols=[0,1,2,3,4],dtype={'names':('filter','lambda','lambda_p','flux','ZP'),'formats':('S10','f8','f8','f8','f8')})

AB_ZPs['lambda']*=1e4 # Transforming wavelength to Angstroms


# -------------SN METADATA ------------------
print " Convert photometry to fluxes to obtain SED "


sn_labels = ['sn','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0']
sn_formats = ['S15','S20','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt('sn_data.dat',dtype={'names':sn_labels,'formats':sn_formats})

current_dir = os.getcwd()
print "Current directory : %s"%current_dir
colormap = pl.cm.spectral



filters = ['W2_uvot','M2_uvot','W1_uvot','U_uvot','B_uvot','V_uvot','U','B','V','R','I','u','g','r','i','z','J','H','K']
filters_colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(filters))]
        
use_filters = ['W2_uvot','M2_uvot','W1_uvot','U','B','V','R','I']

for SN,z_SN,E_B_V,t_0,d_L in SN_DATA[['sn','sn_z','sn_ebv','t_0','hostlumdist']]:
    print SN
    try: 
    
        os.chdir("%s/"%SN)

        # ---------------------LOAD PHOTOMETRY ----------------------------------------------------

        for filter in zip(use_filters):

            filter_file = "mags_%s.out"%filter
            
            if os.path.isfile(filter_file):
                print filter_file,filter
                try:
                    scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = np.genfromtxt(filter_file).T
                    no_nebular = np.where(np.logical_and(scope['jd_%s'%filter]<(t_0+80.),scope['err_%s'%filter]<0.5))[0]
                    scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = \
                                                                                          scope['jd_%s'%filter][no_nebular]-t_0,scope['mag_%s'%filter][no_nebular],scope['err_%s'%filter][no_nebular] 
                                                
                except:
                    
                    print "No errors in band %s ?"%filter
                    scope['jd_%s'%filter],scope['mag_%s'%filter] = np.genfromtxt(filter_file).T
                    no_nebular = np.where(scope['jd_%s'%filter]<(t_0+80.))[0]
                    scope['jd_%s'%filter],scope['mag_%s'%filter] = scope['jd_%s'%filter][no_nebular]-t_0,scope['mag_%s'%filter][no_nebular]
                    
            else:
                print "no %s"%filter_file

        FILTER_ZPS = AB_ZPs[np.in1d(AB_ZPs['filter'],use_filters)]  
        print FILTER_ZPS['filter']
        
        print "# ---------- FILTER_ZPS filters ----------------------------"

        tmin_max,tmax_min = 10000. , -1.0
        max_cadence,min_cadence = 0.,100.0
        min_cadence_filt,max_cadence_filt = None,None
        
        pl.gca().set_color_cycle(filters_colors)
        filters_cadence = []

        for filter in FILTER_ZPS:
            
            
            lco_filter = filter['filter']
            lco_name = lco_filter#lcogt_names[lcogt_filters.index(lco_filter)]
            print lco_filter
            print filter['lambda']
            
            x0 = 1/(filter['lambda']*1e-4)
            R_x = (R_V*a_x(x0)+b_x(x0))
            try:
                scope['mag_%s'%lco_name]
                if len(scope['mag_%s'%lco_name])<=1:
                    continue
            except:
                continue
        
            mag = scope['mag_%s'%lco_name] - R_x*E_B_V   # Deredenned magnitudes 
            err = scope['err_%s'%lco_name]
            t = scope['jd_%s'%lco_name]
            
            pl.errorbar(t,mag,yerr=err,fmt='o',label=lco_filter)

            tmax,tmin = np.max(t),np.min(t)
            dtmax = tmax-tmin
            if tmin >= tmax_min:
                tmax_min = tmin
            if tmax <= tmin_max:
                tmin_max = tmax
        
            scope['min_%s'%lco_filter] ,scope['max_%s'%lco_filter]  = tmin,tmax
            scope['cadence_%s'%lco_filter] = len(t)/dtmax
            scope['inter_%s'%lco_filter] = interp1d(t,mag,kind='linear')
            pl.plot(t,scope['inter_%s'%lco_filter](t),color='k',linestyle='--')
            
            filters_cadence.append(len(t)/dtmax)
            print len(t)/dtmax
            if max_cadence < len(t)/dtmax:
                max_cadence_filt = lco_name
                max_cadence = len(t)/dtmax
        
            if min_cadence > len(t)/dtmax:
                min_cadence_filt = lco_name
                min_cadence = len(t)/dtmax
        
        pl.legend(loc='best',ncol=2,prop={'size':9})
        pl.gca().invert_yaxis()
        pl.show()
    
        baseline = np.arange(int(tmax_min)+1,int(tmin_max),1)
        inter_t = scope['jd_%s'%min_cadence_filt]
        where_inter = np.where(np.logical_and(inter_t <= tmin_max + 0.1 ,inter_t >= tmax_min - 0.1))[0]
        filters_used = []

        for filter in FILTER_ZPS : 
    
            lco_filter = filter['filter']
            lco_name = lco_filter #lcogt_names[lcogt_filters.index(lco_filter)]
            print lco_filter
            x0 = 1/(filter['lambda']*1e-4)
            R_x = (R_V*a_x(x0)+b_x(x0))
    
            try:
                mag = scope['mag_%s'%lco_name] - R_x*E_B_V   # Deredenned magnitudes 
            except:
                continue
                print "Fail"
                
                
            err = scope['err_%s'%lco_name]
            t = scope['jd_%s'%lco_name]
            scope['inter_t_%s'%lco_name] = []
            scope['inter_mags_%s'%lco_name] = []
            scope['inter_err_%s'%lco_name] = []
            filters_used.append(lco_name)
            
            for t_inter in inter_t[where_inter]:
                
                nearest_t = np.where(np.abs(t-t_inter)<0.25)[0]
         
                if len(nearest_t)>0:

                    scope['inter_t_%s'%lco_name].append(t[nearest_t[0]])
                    scope['inter_mags_%s'%lco_name].append(mag[nearest_t[0]])
                    scope['inter_err_%s'%lco_name].append(err[nearest_t[0]])
            
                else:
                    
                    scope['inter_t_%s'%lco_name].append(t_inter)
                    inter,inter_err = inter_sample(t_inter,t,mag,yerr=err,n_sample=1000,plot=False,confidence=False)
                    scope['inter_mags_%s'%lco_name].append(inter)
                    scope['inter_err_%s'%lco_name].append(inter_err)
                    
                    # interpolar el resto de los filtros en la misma epoca    

            
        inter_mags = [ (filter,np.array(scope['inter_mags_%s'%filter]),np.array(scope['inter_err_%s'%filter])) \
                   for filter in filters_used]

        sorted_mags = sorted(inter_mags,\
                                key = lambda x : FILTER_ZPS['lambda'][np.where(FILTER_ZPS['filter']==x[0])[0]])

        sorted_lambdas = np.sort(FILTER_ZPS['lambda_p'][np.in1d(FILTER_ZPS['filter'],filters_used)])
        sorted_bands = [FILTER_ZPS['filter'][np.where(FILTER_ZPS['lambda_p']==lam)[0]][0] for lam in sorted_lambdas]

        fluxes = [[] for i in range(len(where_inter))]
        fluxes_err = [[] for i in range(len(where_inter))]
        
        for filter,mags,errs in sorted_mags:
            
            print filter
            index = filters_used.index(filter)
            lambda_f = FILTER_ZPS['lambda_p'][index]
            ZP = FILTER_ZPS['ZP'][index]
            print lambda_f,ZP
            
            # m_AB = -2.5log(f_AB) - 48.6 - ZP
            F_AB = 10**(-0.4*(48.6+mags+ZP)) # Jansky 
            F_err =  F_AB*0.4*errs*np.log(10)
            
            pl.errorbar(inter_t[where_inter],mags,yerr = errs,label=filter)
    
            for i,f in enumerate(zip(F_AB,F_err)):
                
                fab = f[0]
                ferr = f[1]
                fluxes[i].append(fab)
                fluxes_err[i].append(np.asarray(ferr))
                
        
        pl.legend()
        pl.gca().invert_yaxis()
        pl.show()
        
        
        print "# --------- Saving SEDs -------------------------------------"

        d_SED = {}
        d_SED['t'] = inter_t[where_inter]
        d_SED['data'] = fluxes
        d_SED['daterr'] = fluxes_err
        d_SED['lambda'] = sorted_lambdas
        d_SED['bands'] = sorted_bands 
        np.savez( 'SED_%s.npz'%SN, **d_SED )
        
        
        print "# --------- Plot SED's --------------------------------------"
    

        fig = pl.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        
        colormap = pl.cm.spectral
        fig.gca().set_color_cycle([colormap(i) for i in np.linspace(0.1, 0.9,len(where_inter))])
        
        filter_names = AB_ZPs['filter'][np.in1d(sorted_lambdas,AB_ZPs['lambda'])]
        responses =  lcogt_responses
        
        for epoch,flux,flux_err in zip(inter_t[where_inter],fluxes,fluxes_err):
            
            flux = np.asarray(flux)
            ax.errorbar(c_cms/(sorted_lambdas*1e-8),flux,yerr=flux_err,linestyle='--',label='%2.2f'%epoch)
                
                
        ax.set_title(SN)
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('Flux in ergs/s/cm2/Hz')
        ax.legend(loc='best',ncol=2,prop={'size':9})
        fig.savefig('./SED_nu_%s.png'%R_V)
        pl.show()
    
        fig = pl.figure(figsize=(10,6))
        ax = fig.add_subplot(111)

        colormap = pl.cm.spectral
        fig.gca().set_color_cycle([colormap(i) for i in np.linspace(0.1, 0.9, len(where_inter))])

        for epoch,flux,flux_err in zip(inter_t[where_inter],fluxes,fluxes_err):
            
            flux = np.asarray(flux)
            ax.errorbar(sorted_lambdas,flux*c_cms/((sorted_lambdas*1e-4)**2.),\
                        yerr=c_cms*np.asarray(flux_err)/((sorted_lambdas*1e-4)**2.),marker='o',linestyle='--',label='%2.2f'%epoch)
            
        trans = ax.get_xaxis_transform()
        for lam,flux in zip(sorted_lambdas,fluxes[0]):
            flux = np.asarray(flux)
            max_flux = np.max(flux)*c_cms/((lam*1e-4)**2.)
            band = FILTER_ZPS['filter'][np.where(FILTER_ZPS['lambda_p']==lam)[0]]
            print band,lam,max_flux
            ax.axvline(lam,linestyle='--',color='k',alpha=0.4)
            ax.annotate(str(band[0][0:2]),(lam+5,1.05),xycoords=trans,size=8,rotation=90)
                
        ax.set_title(SN)
        ax.set_xlabel(r'$\mathrm{Rest \ Wavelength \ [\AA]}$')
        ax.set_ylabel('Flux in ergs/s/cm2/A')
        ax.legend(loc='best',ncol=2,prop={'size':9})
        fig.savefig('./SED_lam_%s.png'%(str(R_V)))
        pl.show()
        pl.close()
        
        print "------------ Integrating SED -------------------------"
        
        fig = pl.figure(figsize=(10,6))
        ax = fig.add_subplot(111)

        from scipy.integrate import simps
        L_bol = []
        M_bol = []
        Lfile = open('%s_Lbol.dat'%SN,'w')
        
        for epoch,flux,flux_err in zip(inter_t[where_inter],fluxes,fluxes_err):
            
            y,x = flux,c_cms/(sorted_lambdas*1e-8) 
            I = -simps(y,x)
            L = 4*pi * I * (d_L*1e6*pc_to_cm)**2.
            print epoch,I,4*pi*I*(d_L*1e6*pc_to_cm)**2. 
            L_bol.append(L)
            ax.plot(epoch,L,marker='o',color='k')
            Lfile.write("%s \t %s \n"%(epoch,L))
            
        Lfile.close()
        ax.set_xlabel('t')
        ax.set_ylabel('L  [erg/s]')
        fig.savefig('./Lbol_%s.png'%(SN))
        pl.show()
        pl.close()
        
        
        try:
            os.chdir(current_dir)
            
        except:
            
            print "No diretory %s"%current_dir

    except:
        
        print "No directory %s/"%SN
        continue


        
try: 

    os.chdir(current_dir)

except:

    print "No diretory %s"%current_dir
