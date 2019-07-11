from astro_utils import *

scope = globals()

# LCOGT 

LCOGT_filters = ('bssl-bx-004','bssl-vx-022','bssl-rx-007','bssl-ix-023','SDSS.gp','SDSS.rp','SDSS.ip')
LCOGT_filters_id = ('lco_B','lco_V','lco_R','lco_I')
LCOGT_colors = ('b','g','y','r')

LCOGT_DIR = "/home/dust_speck/SN/lcogt_filters"
f_name = ('lambda','response')
f_format = ('f8','f8')

for filter_name,filter_id in zip(LCOGT_filters,LCOGT_filters_id):

    scope['%s_response'%filter_id] = np.genfromtxt("%s/%s.txt"%(LCOGT_DIR,filter_name),dtype={'names':f_name,'formats':f_format})
    scope['%s_response'%filter_id]['lambda'] *=10
# SWIFT 

UVOT = ('UVW2','UVM2','UVW1','U','B','V')
UVOT_names = ('UVW2','UVM2','UVW1','U')
UVOT_colors = ("indigo",'blueviolet','mediumpurple','darkslateblue')
UVOT_DIR = "/home/dust_speck/SN/SWIFT/filters"
f_name = ('lambda','response')
f_format = ('f8','f8')

for uvot_filter,uvot_name in zip(UVOT,UVOT_names):

    scope['%s_response'%uvot_name] = np.genfromtxt("%s/Swift_UVOT.%s.dat"%(UVOT_DIR,uvot_filter),dtype={'names':f_name,'formats':f_format})
    scope['%s_response'%uvot_name]['response']*=10/659.
# ----2MASS-------------------------------------------

twoMASS = ('2MASS.J','2MASS.H','2MASS.Ks')
twoMASS_names = ('J','H','Ks')
twoMASS_colors = ('firebrick','maroon','darkred')
twoMASS_DIR = "/home/dust_speck/SN/2MASS_filters"
f_name = ('lambda','response')
f_format = ('f8','f8')

for uvot_filter,uvot_name in zip(twoMASS,twoMASS_names):

    scope['%s_response'%uvot_name] = np.genfromtxt("%s/%s.dat"%(twoMASS_DIR,uvot_filter),dtype={'names':f_name,'formats':f_format})


# ---------------------------------------------

pl.figure(figsize=(20,10))

T_bb = 5000 # BB temperature in kelvins

lamb = np.arange(1000,30000,1)
B_lamb =  B_lambda(lamb*1e-8,T=T_bb)*1e-8 ## erg/s/cm2/Ang
pl.plot(lamb,B_lamb,color='k')

der_B_lamb = lambda y : derivative(lambda x: B_lambda(x*1e-8,T=T_bb)*1e-8,y,n=1,dx = 0.1)
der2_B_lamb = lambda y : derivative(lambda x: B_lambda(x*1e-8,T=T_bb)*1e-8,y,n=2,dx = 0.1)

B_max = np.max(B_lamb)

Data = {'lambda':lamb,'Flux':B_lamb}

offsets= [[],[]]
filter_lambdas = []

for filter_id,color in zip(UVOT_names+LCOGT_filters_id+twoMASS_names,UVOT_colors+LCOGT_colors+twoMASS_colors):
    
    lambdas,response = scope["%s_response"%filter_id]['lambda'],scope["%s_response"%filter_id]['response']
    scope["lambda_cen_%s"%filter_id]  = simps(response*lambdas,lambdas)/simps(response,lambdas)
    lam_c = scope["lambda_cen_%s"%filter_id]
    scope["width_%s"%filter_id]  = simps(response*(lambdas-lam_c)**2.,lambdas)/simps(response,lambdas)
    width_lam = scope["width_%s"%filter_id]
    pl.axvline(lam_c,color=color,linestyle='-')
    pl.fill_between(lambdas,0,response*B_max,color=color,alpha=0.4)
    pl.annotate(filter_id,xy=(lam_c+50,0.4*B_max),xycoords='data',color=color,rotation=90,size=15)    
    F_s = synthetic_mag(Data,scope["%s_response"%filter_id],lambda_=True,mag=False,ZP=0.0,F_0=363.1e-11,z=0.0)
    F_c = B_lambda(lam_c*1e-8,T=T_bb)*1e-8
    error_F = (F_s - F_c)/width_lam
    offsets[0].append(lam_c)
    offsets[1].append(error_F)
    pl.plot(lam_c,F_s,marker='o',color=color,markersize=10)
    filter_lambdas.append(lam_c)
    print filter_id,lam_c,width_lam,F_s,error_F 


pl.ylabel(r'$F_{\lambda} \ \mathrm{[erg/s/cm^2/\AA]}$',size=20)
pl.xlabel(r'$ \lambda \ \mathrm{[\AA]}$',size=20)
pl.annotate(r'$T = %4.1f \ [K]$'%T_bb,xy=(9200, 0.9*B_max),xycoords='data',size=20)
pl.xticks(np.arange(min(lamb), max(lamb)+1, 1000))
pl.tight_layout()
pl.xlim(1000,)
pl.savefig('/home/dust_speck/Descargas/BB_plot.png')
pl.show()

pl.plot(offsets[0],np.asarray(offsets[1]),marker='o',label=r'$\left(F_s - F_{\lambda_0}\right)/\left<\Delta\lambda^2\right>$')
pl.plot(lamb,der2_B_lamb(lamb)/2.,label=r'$\frac{1}{2}\frac{d^2F_{\lambda_0}}{d\lambda^2}$')
pl.xticks(np.arange(min(lamb), max(lamb)+1, 1000))
pl.tight_layout()
pl.xlabel(r'$ \lambda \ \mathrm{[\AA]}$',size=20)
pl.legend(prop={'size':20})
pl.show()

T_dust = 3500
R_V = 3.1
E_B_Vs = [0.1,0.25,0.5]

fig1 = pl.figure()
ax = fig1.add_subplot(111)

R_lam = lambda x: R_V*a_x(x)+b_x(x)

for E_B_V in E_B_Vs:
    
    x_lam = 1/(lamb*1e-4)
    filter_lambdas = np.asarray(filter_lambdas)
    tau_lam = np.log(10)*E_B_V*np.array([R_lam(x_) for x_ in x_lam])/2.5
    print tau_lam
    I_lam = B_lambda(lamb*1e-8,T=T_bb)*np.exp(-tau_lam) + (1-np.exp(-tau_lam))*B_lambda(lamb*1e-8,T=T_dust)
    T_lam = T_BB(I_lam,lamb*1e-8)
    plot = ax.plot(lamb,T_lam)
    color = plot[0].get_color()
    
    x_filt = 1/(filter_lambdas*1e-4)
    tau_filt = np.log(10)*E_B_V*np.array([R_lam(x_) for x_ in x_filt])/2.5
    I_filt = B_lambda(filter_lambdas*1e-8,T=T_bb)*np.exp(-tau_filt) + (1-np.exp(-tau_filt))*B_lambda(filter_lambdas*1e-8,T=T_dust)
    T_filt = T_BB(I_filt,filter_lambdas*1e-8)
    ax.plot(filter_lambdas,T_filt,marker='o',linestyle='None',label=r'$E(B-V) = %2.2f$'%E_B_V,color=color)
    ax.plot(lamb,T_bb*np.exp(-tau_lam)+T_dust*(1-np.exp(-tau_lam)),linestyle='--',color=color)
    


ax.axhline(T_bb,color='k')
ax.axhline(T_dust,color='k')    
pl.legend(loc='best')
pl.xlabel(r"$\lambda$",size=15)
pl.ylabel(r"$T_B$",size=15)
pl.xlim(1200,)
pl.show()

pl.plot(filter_lambdas,I_filt)
pl.plot(lamb,I_lam,label="I_lam(%s)"%E_B_V)
pl.plot(lamb,B_lamb*1e8,label='T_*')
pl.plot(lamb,B_lambda(lamb*1e-8,T=T_dust),label='T_n')
pl.legend()
pl.show()


pl.plot(lamb,T_bb*np.exp(-tau_lam)+T_dust*(1-np.exp(-tau_lam)))
