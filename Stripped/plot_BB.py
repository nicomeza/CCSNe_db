from astro_utils import *

scope = globals()

# LCOGT 

LCOGT_filters = ('bssl-bx-004','bssl-vx-022','bssl-rx-007','bssl-ix-023','SDSS.gp','SDSS.rp','SDSS.ip')
LCOGT_filters_id = ('lco_B','lco_V','lco_R','lco_I','SDSS.gp','SDSS.rp','SDSS.ip')
LCOGT_colors = ('b','g','y','r')

LCOGT_DIR = "/home/dust_speck/SN/lcogt_filters"
f_name = ('lambda','response')
f_format = ('f8','f8')

for filter_name,filter_id in zip(LCOGT_filters,LCOGT_filters_id):

    scope['%s_response'%filter_id] = np.genfromtxt("%s/%s.txt"%(LCOGT_DIR,filter_name),dtype={'names':f_name,'formats':f_format})
    s_x,s_y = (10. * scope['%s_response'%filter_id]['lambda'],scope['%s_response'%filter_id]['response'])
    where = np.where(s_y>0)[0]

# SWIFT 

UVOT = ('UVW2','UVM2','UVW1','U','B','V')
UVOT_names = ('UVW2','UVM2','UVW1','U','B','V')
UVOT_colors = ("indigo",'blueviolet','mediumpurple','darkslateblue')
UVOT_DIR = "/home/dust_speck/SN/SWIFT/filters"
f_name = ('lambda','response')
f_format = ('f8','f8')

for uvot_filter,uvot_name in zip(UVOT,UVOT_names):

    scope['%s_response'%uvot_name] = np.genfromtxt("%s/Swift_UVOT.%s.dat"%(UVOT_DIR,uvot_filter),dtype={'names':f_name,'formats':f_format})
    s_x,s_y = (scope['%s_response'%uvot_name]['lambda'],scope['%s_response'%uvot_name]['response'])
    where = np.where(s_y>0)[0]
    s_x,s_y = s_x[where],s_y[where]

# ----2MASS-------------------------------------------

twoMASS = ('2MASS.J','2MASS.H','2MASS.Ks')
twoMASS_names = ('J','H','Ks')
twoMASS_colors = ('firebrick','maroon','darkred')
twoMASS_DIR = "/home/dust_speck/SN/2MASS_filters"
f_name = ('lambda','response')
f_format = ('f8','f8')

for uvot_filter,uvot_name in zip(twoMASS,twoMASS_names):

    scope['%s_response'%uvot_name] = np.genfromtxt("%s/%s.dat"%(twoMASS_DIR,uvot_filter),dtype={'names':f_name,'formats':f_format})
    s_x,s_y = (scope['%s_response'%uvot_name]['lambda'],scope['%s_response'%uvot_name]['response'])
    where = np.where(s_y>0)[0]
    s_x,s_y = s_x[where],s_y[where]


# ---------------------------------------------

pl.figure(figsize=(20,10))

T_BB = 5000 # BB temperature in kelvins

lamb = np.arange(1000,25000,1)
B_lamb =  B_lambda(lamb*1e-8,T=T_BB)*1e-8 ## erg/s/cm2/Ang
pl.plot(lamb,B_lamb,color='k')

B_max = np.max(B_lamb)

for filter_id,color in zip(LCOGT_filters_id,LCOGT_colors):
    
    lambdas,response = scope["%s_response"%filter_id]['lambda'],scope["%s_response"%filter_id]['response']
    scope["lambda_cen_%s"%filter_id]  = simps(response*lambdas,lambdas)/simps(response,lambdas)
    pl.axvline(10*scope["lambda_cen_%s"%filter_id],color=color,linestyle='-')
    pl.fill_between(lambdas*10,0,response*B_max,color=color,alpha=0.4)
    pl.annotate(filter_id,xy=(10*scope["lambda_cen_%s"%filter_id]+50,0.4*B_max),xycoords='data',color=color,rotation=90,size=15)    
    print filter_id,scope["lambda_cen_%s"%filter_id] 


for filter_id,color in zip(UVOT_names,UVOT_colors):
    
    lambdas,response = scope["%s_response"%filter_id]['lambda'],scope["%s_response"%filter_id]['response']
    scope["lambda_cen_%s"%filter_id]  = simps(response*lambdas,lambdas)/simps(response,lambdas)
    pl.axvline(scope["lambda_cen_%s"%filter_id],color=color,linestyle='-')
    pl.fill_between(lambdas,0,response*B_max*10/659.,color=color,alpha=0.4)
    pl.annotate(filter_id,xy=(scope["lambda_cen_%s"%filter_id]+50,0.4*B_max),xycoords='data',color=color,rotation=90,size=15)    
    print filter_id,scope["lambda_cen_%s"%filter_id] 


for filter_id,color in zip(twoMASS_names,twoMASS_colors):
    
    lambdas,response = scope["%s_response"%filter_id]['lambda'],scope["%s_response"%filter_id]['response']
    scope["lambda_cen_%s"%filter_id]  = simps(response*lambdas,lambdas)/simps(response,lambdas)
    pl.axvline(scope["lambda_cen_%s"%filter_id],color=color,linestyle='-')
    pl.fill_between(lambdas,0,response*B_max,color=color,alpha=0.4)
    pl.annotate(filter_id,xy=(scope["lambda_cen_%s"%filter_id]+50,0.4*B_max),xycoords='data',color=color,rotation=90,size=15)    
    print filter_id,scope["lambda_cen_%s"%filter_id] 


pl.ylabel(r'$F_{\lambda} \ \mathrm{[erg/s/cm^2/\AA]}$',size=20)
pl.xlabel(r'$ \lambda \ \mathrm{[\AA]}$',size=20)
pl.annotate(r'$T = %4.3f \ [K]$'%T_BB,xy=(9200, 0.9*B_max),xycoords='data',size=20)
pl.xticks(np.arange(min(lamb), max(lamb)+1, 1000))
pl.tight_layout()
pl.xlim(1000,)
pl.savefig('/home/dust_speck/Descargas/BB_plot.png')
pl.show()



