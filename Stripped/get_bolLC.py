import numpy as np
import matplotlib.pyplot as pl
import astro_utils


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



# Convert to fluxes to obtain SED 

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

UVOT = ('uvw2','uvm2','uvw1','U','B','V')
UVOT_names = ('W2','M2','W1','U','B','V')
UVOT_DIR = "/home/dust_speck/SN/SWIFT"
f_name = ('lambda','response')
f_format = ('f8','f8')
uvot_responses = []

for uvot,uvot_name in zip(UVOT,UVOT_names):
    
    scope['%s_response'%uvot_name] = np.genfromtxt("%s/%s_swift.txt"%(UVOT_DIR,uvot),dtype={'names':f_name,'formats':f_format})
    s_x,s_y = (scope['%s_response'%uvot_name]['lambda'],scope['%s_response'%uvot_name]['response'])
    where = np.where(s_y>0)[0]
    uvot_responses.append((s_x[where]*1e-8,s_y[where]))

# ---------------------  SN METADATA ---------------------------------------------------

z_14jb = 0.006031
SN = "ASASSN-14jb"
E_B_V_14jb = 0.0154 # Local Reddening

# ---------------------LOAD PHOTOMETRY ---------------------------------------------------




filters = {'UWW2','UWM2','UWW1','U','B','V'}
lcogt_filters = ['B','V','g','r','i']

for filter in zip(filters):

    filter_file = "mags_%s.out"%filter
    print filter_file,filter
    scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = np.genfromtxt(filter_file).T
    no_nebular = np.where(scope['jd_%s'%filter]<80.)[0]
    scope['jd_%s'%filter],scope['mag_%s'%filter],scope['err_%s'%filter] = \
        scope['jd_%s'%filter][no_nebular],scope['mag_%s'%filter][no_nebular],scope['err_%s'%filter][no_nebular] 
    #print scope['jd_%s'%filter][-1]

LCOGT = LCOGT_data[np.in1d(LCOGT_data['filter'],lcogt_filters[2:])]  

print LCOGT['filter']
    

# --------Optical Spectra ----------------------------------------------------------

read_spectra = False
if read_spectra:
    
    list = "../spectra/spectra.dat"
    spectras = np.genfromtxt(list,usecols=[0,1],dtype={'names':('fits','JD'),'formats':('S35','f8')})
    spectras['JD'] -= (t_0-2400000.5)


#------------------ UVOT ZERO POINT -------------------------------------

uvot = np.genfromtxt('/home/dust_speck/SN/UVOT_zeropoint.dat',usecols=[0,1,2],dtype={'names':('filter','lambda','flux'),'formats':('S5','f8','f8')})
vega = np.genfromtxt('/home/dust_speck/SN/Vega_zeropoints.dat',usecols=[0,1,2],dtype={'names':('filter','lambda','flux'),'formats':('S5','f8','f8')})
uvot_AB = np.genfromtxt('/home/dust_speck/SN/UVOT_AB.dat',usecols=[0,1,2,3],dtype={'names':('filter','z','z_err','AB_V'),'formats':('S5','f8','f8','f8')})
#print uvot_AB['filter'],uvot_AB['AB_V']
tmin_max,tmax_min = 10000. , -1.0
tmin_min,tmax_max = 19999, -0.1
max_cadence,min_cadence = 0.,100.0
min_cadence_filt,max_cadence_filt = None,None
tmin_filt = None

print "# --------------------uvot flux conversion ---------------------------------"


for filter in uvot:
    
    uvot_filter = filter['filter']
    filter_name = uvot_names[uvot_filters.index(uvot_filter)]
    #print filter_name,uvot_filter
    #print filter['lambda']
    x0 = 1/(filter['lambda']*1e-4)
    R_x = (R_V*a_x(x0)+b_x(x0))
    #print "R_%s = %s"%(filter['filter'],R_x)
    mag = scope['mag_%s'%filter_name] - R_x*E_B_V_14jb   # Deredenned magnitudes 
    err = scope['err_%s'%filter_name]
    t = scope['jd_%s'%filter_name] 
    
    pl.errorbar(t,mag,yerr=err,fmt='o',label=uvot_filter)

    tmax,tmin = np.max(t),np.min(t)
    dtmax = tmax-tmin
    
    if tmin >= tmax_min:
        tmax_min = tmin
    if tmin<=tmin_min:
        tmin_min = tmin
        tmin_filt  = filter_name
    if tmax <= tmin_max:
        tmin_max = tmax
    if tmax > tmax_max:
        tmax_max = tmax
        
    scope['min_%s'%filter_name] ,scope['max_%s'%filter_name]  = tmin,tmax
    scope['cadence_%s'%filter_name] = len(t)/dtmax
    scope['inter_%s'%filter_name] = interp1d(t,mag,kind='linear')
    pl.plot(t,scope['inter_%s'%filter_name](t),label='%s inter'%uvot_filter)
    
    if max_cadence < len(t)/dtmax:
        max_cadence_filt = filter_name
        max_cadence = len(t)/dtmax
        
    if min_cadence > len(t)/dtmax:
        min_cadence_filt = filter_name
        min_cadence = len(t)/dtmax
        
pl.legend()
pl.gca().invert_yaxis()
pl.show()
    
print tmin_min,tmax_max,tmin_filt
baseline = np.arange(int(tmax_min)+1,int(tmin_max),1)
inter_t = scope['jd_%s'%tmin_filt]
#print inter_t
where_inter = np.where(np.logical_and(inter_t <= tmax_max + 0.1 ,inter_t>= tmin_min - 0.1))[0]
#print min_cadence_filt,min_cadence
#print max_cadence_filt, max_cadence


for filter in uvot : 
    
    uvot_filter = filter['filter']
    filter_name = uvot_names[uvot_filters.index(uvot_filter)]
    x0 = 1/(filter['lambda']*1e-4)
    R_x = (R_V*a_x(x0)+b_x(x0))
    #print "R_%s = %s"%(filter['filter'],R_x)
    mag = scope['mag_%s'%filter_name] - R_x*E_B_V_14jb   # Deredenned magnitudes 
    err = scope['err_%s'%filter_name]
    print filter_name, t 
    t = scope['jd_%s'%filter_name] 
    scope['inter_t_%s'%filter_name] = []
    scope['inter_mags_%s'%filter_name] = []
    scope['inter_err_%s'%filter_name] = []
    
    for t_inter in inter_t[where_inter]:
        
        try:
            
            nearest_t = np.where(np.abs(t-t_inter)<0.1)[0]
            print np.where(np.abs(t-t_inter)<0.1)
            if len(nearest_t)>0:
                nearest_t = nearest_t[0]
                scope['inter_t_%s'%filter_name].append(t[nearest_t])
                scope['inter_mags_%s'%filter_name].append(mag[nearest_t])
                scope['inter_err_%s'%filter_name].append(err[nearest_t])
            
            else:
                print "interpolated mag in %s band at %s +- 0.1"%(filter_name,t_inter)
                scope['inter_t_%s'%filter_name].append(t_inter)
                inter,inter_err = inter_sample(t_inter,t,mag,yerr=err,n_sample = 1000,plot=False)
                scope['inter_mags_%s'%filter_name].append(inter)
                scope['inter_err_%s'%filter_name].append(inter_err)
        except:
            
            print "no mag in %s band at %s +- 0.1"%(filter_name,t_inter)
            scope['inter_t_%s'%filter_name].append(t_inter)
            scope['inter_mags_%s'%filter_name].append(99)
            scope['inter_err_%s'%filter_name].append(99)
            

                
# interpolar el resto de los filtros en la misma epoca    

inter_mags = [ (filter,np.array(scope['inter_mags_%s'%filter]),np.array(scope['inter_err_%s'%filter])) \
                  for filter in  uvot_names if filter != tmin_filt]

inter_mags.append((tmin_filt,scope['mag_%s'%tmin_filt][where_inter],scope['err_%s'%tmin_filt][where_inter]))
sorted_mags = sorted(inter_mags,\
         key = lambda x : uvot['lambda'][np.where(uvot['filter']==uvot_filters[uvot_names.index(x[0])])[0]])

sorted_lambdas_uvot = np.asarray(sorted(uvot['lambda']))
print sorted_mags,sorted_lambdas_uvot

fluxes_uvot = [[] for i in range(len(where_inter))]
fluxes_uvot_err = [[] for i in range(len(where_inter))]

for filter,mags,errs in sorted_mags:
    
    print filter
    index = np.where(uvot['filter']==uvot_filters[uvot_names.index(filter)])[0]
    #print uvot['filter'][index]
    lambda_f = uvot['lambda'][index]
    #print lambda_f
    z,zerr = uvot_AB['z'][index],uvot_AB['z_err'][index]
    dab = uvot_AB['AB_V'][index]
    #print dab
    F_AB = 10**(-0.4*(48.6+mags+dab)) # Jansky
    F_err = F_AB*0.4*errs*np.log(10) 

    
    pl.errorbar(inter_t[where_inter],mags,yerr = errs,label=filter)
    
    for i,f in enumerate(zip(F_AB,F_err)):
        
        fab = f[0]
        ferr = f[1]
        fluxes_uvot[i].append(fab)
        fluxes_uvot_err[i].append(np.asarray(ferr))


pl.legend()
pl.gca().invert_yaxis()
pl.show()


tmin_max,tmax_min = 10000. , -1.0
max_cadence,min_cadence = 0.,100.0
min_cadence_filt,max_cadence_filt = None,None

    
print "# ---------- LCOGT filters ----------------------------"


for filter in LCOGT:
    
    lco_filter = filter['filter']
   
    print lco_filter
    print filter['lambda']
    
    x0 = 1/(filter['lambda']*1e-4)
    R_x = (R_V*a_x(x0)+b_x(x0))
    #print "R_%s = %s"%(filter['filter'],R_x)
    mag = scope['mag_%s'%lco_filter] - R_x*E_B_V_14jb   # Deredenned magnitudes 
    err = scope['err_%s'%lco_filter]
    t = scope['jd_%s'%lco_filter] 
    
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
    pl.plot(t,scope['inter_%s'%lco_filter](t),label='%s inter'%lco_filter)
    
    if max_cadence < len(t)/dtmax:
        max_cadence_filt = lco_filter
        max_cadence = len(t)/dtmax
        
    if min_cadence > len(t)/dtmax:
        min_cadence_filt = lco_filter
        min_cadence = len(t)/dtmax
        
pl.legend()
pl.gca().invert_yaxis()
pl.show()
    

baseline = np.arange(int(tmax_min)+1,int(tmin_max),1)
inter_t_lco = scope['jd_%s'%min_cadence_filt]
where_inter_lco = np.where(np.logical_and(inter_t_lco <= tmin_max + 0.1 ,inter_t_lco >= tmax_min - 0.1))[0]
#print min_cadence_filt,min_cadence
#print max_cadence_filt, max_cadence


for filter in LCOGT : 
    
    lco_filter = filter['filter']
    print lco_filter
    x0 = 1/(filter['lambda']*1e-4)
    R_x = (R_V*a_x(x0)+b_x(x0))
    
    mag = scope['mag_%s'%lco_filter] - R_x*E_B_V_14jb   # Deredenned magnitudes 
    err = scope['err_%s'%lco_filter]
    t = scope['jd_%s'%lco_filter] 
    scope['inter_t_%s'%lco_filter] = []
    scope['inter_mags_%s'%lco_filter] = []
    scope['inter_err_%s'%lco_filter] = []
    
    for t_inter in inter_t[where_inter]:
        
        nearest_t = np.where(np.abs(t-t_inter)<0.25)[0]
         
        if nearest_t:
            
            #print  t[nearest_t[0]],mag[nearest_t[0]],err[nearest_t[0]]
            scope['inter_t_%s'%lco_filter].append(t[nearest_t[0]])
            scope['inter_mags_%s'%lco_filter].append(mag[nearest_t[0]])
            scope['inter_err_%s'%lco_filter].append(err[nearest_t[0]])
            
        else:
            
            scope['inter_t_%s'%lco_filter].append(t_inter)
            inter,inter_err = inter_sample(t_inter,t,mag,yerr=err,n_sample=1000,plot=False)
            scope['inter_mags_%s'%lco_filter].append(inter)
            scope['inter_err_%s'%lco_filter].append(inter_err)
            #print t_inter,inter,inter_err
            
    


                
# interpolar el resto de los filtros en la misma epoca    

#inter_mags_lco = [ (filter,np.array(scope['inter_mags_%s'%filter]),np.array(scope['inter_err_%s'%filter])) \
#                  for filter in  LCOGT['filter'] if filter != min_cadence_filt]

inter_mags_lco = [ (filter,np.array(scope['inter_mags_%s'%filter]),np.array(scope['inter_err_%s'%filter])) \
                  for filter in  LCOGT['filter'] ]

#inter_mags_lco.append((min_cadence_filt,scope['mag_%s'%min_cadence_filt][where_inter_lco],2.*scope['err_%s'%min_cadence_filt][where_inter_lco]))

sorted_mags_lco = sorted(inter_mags_lco,\
         key = lambda x : LCOGT['lambda'][np.where(LCOGT['filter']==x[0])[0]])

sorted_lambdas_lco = np.sort(LCOGT['lambda'])
#print sorted_mags_lco,sorted_lambdas_lco

fluxes_lco = [[] for i in range(len(where_inter))]
fluxes_lco_err = [[] for i in range(len(where_inter))]

for filter,mags,errs in sorted_mags_lco:
    
    print filter
    index = np.where(LCOGT['filter']==filter)[0]
    print LCOGT['filter'][index]
    lambda_f = LCOGT['lambda'][index]
    print lambda_f
    
    F_AB = 10**(-0.4*(48.6+mags)) # Jansky
    F_err =  F_AB*0.4*errs*np.log(10)
    
    
    pl.errorbar(inter_t[where_inter],mags,yerr = errs,label=filter)
    
    for i,f in enumerate(zip(F_AB,F_err)):
        
        fab = f[0]
        ferr = f[1]
        fluxes_lco[i].append(fab)
        fluxes_lco_err[i].append(np.asarray(ferr))


pl.legend()
pl.gca().invert_yaxis()
pl.show()


print "# --------- Plot SED's --------------------------------------"
    

fig = pl.figure(figsize=(10,6))
ax = fig.add_subplot(111)

colormap = pl.cm.spectral
fig.gca().set_color_cycle([colormap(i) for i in np.linspace(0.1, 0.9,len(where_inter))])

fluxes = np.concatenate([fluxes_uvot,fluxes_lco],axis=1)
sorted_lambdas = np.concatenate([sorted_lambdas_uvot,sorted_lambdas_lco],axis=0)
fluxes_err = np.concatenate([fluxes_uvot_err,fluxes_lco_err],axis=1)
responses =  np.concatenate([uvot_responses,lcogt_responses])
#print all_fluxes,all_lambdas


for epoch,flux,flux_err in zip(inter_t[where_inter],fluxes,fluxes_err):
    
    flux = np.asarray(flux)
    ax.errorbar(c_cms/(sorted_lambdas*1e-8),flux,yerr=flux_err,linestyle='--',label='%2.2f'%epoch)
    

    
ax.set_xlabel('Frequency [Hz]')
ax.set_ylabel('Flux in ergs/s/cm2/Hz')
#ax.set_xlim(0.2,1.9)
ax.legend(loc='best',ncol=2,prop={'size':9})
fig.savefig('/home/dust_speck/SED_nu_%s.png'%R_V)
pl.show()
    
fig = pl.figure(figsize=(10,6))
ax = fig.add_subplot(111)

colormap = pl.cm.spectral
fig.gca().set_color_cycle([colormap(i) for i in np.linspace(0.1, 0.9, len(fluxes)+len(fluxes_lco))])




for epoch,flux,flux_err in zip(inter_t[where_inter],fluxes,flux_err):
    
    print flux_err/flux
    flux = np.asarray(flux)
    ax.errorbar(sorted_lambdas,flux*c_cms/((sorted_lambdas*1e-4)**2),\
                yerr=c_cms*np.asarray(flux_err)/((sorted_lambdas*1e-4)**2),linestyle='--',label='%2.2f'%epoch)
    




ax.set_xlabel(r'$\mathrm{Rest \ Wavelength \ [\AA]}$')
ax.set_ylabel('Flux in ergs/s/cm2/A')
ax.legend(loc='best',ncol=2,prop={'size':9})
fig.savefig('/home/dust_speck/SED_lam_%s.png'%(R_V))
pl.show()
pl.close()
    