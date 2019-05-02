#-*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.ticker as ticker
import os 
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d,splrep,splev,UnivariateSpline
import sgolay
from scipy.integrate import simps
from scipy.special import erf
import scipy.integrate as integrate
from matplotlib.patches import Ellipse

def overlap(a, b):
    # return the indices in a that overlap with b, also returns 
    # the corresponding index in b only works if both a and b are unique! 
    # This is not very efficient but it works
    bool_a = np.in1d(a,b)
    ind_a = np.arange(len(a))
    ind_a = ind_a[bool_a]

    ind_b = np.array([np.argwhere(b == a[x]) for x in ind_a]).flatten()
    return ind_a,ind_b


def linear_model(x,a,b):
    return a*x+b

def bol_correction(color,color_mag): # From Pejcha-Prieto 2015
    
    if color == "B-V":  # valid range:(−0.26, 1.27). rms = 0.023 
        
        c_0 = -0.3716
        c_1 = 2.9669
        c_2 = -6.2797
        c_3 = 5.8950
        c_4 = -2.0233
        #rms = 0.023
        
        BC = np.polynomial.polynomial.polyval(color_mag,[c_0,c_1,c_2,c_3,c_4])
        
        if np.logical_or(color_mag > 1.27,color_mag < -0.26).any():
            print "Warning : Out of range Color found"
        
        return BC
        
    else:
        print "The first input is a string with the color used."
        return 0.

def bol_correction_2(color,color_mag): # From Lyman 2014
    
    if color == "B-V":  # valid range:(0.0, 1.6). rms = 0.094
        
        c_0 = -0.138
        c_1 = -0.013
        c_2 = -0.649

        BC = np.polynomial.polynomial.polyval(color_mag,[c_0,c_1,c_2])
              
        if np.logical_or(color_mag > 1.6,color_mag < 0.0).any():
            print "Warning : Out of range Color found"
        
        return BC
    
    if color == "g-i":  # valid range:(-0.5, 1.4). rms = 0.022
        
        c_0 = -0.007
        c_1 = -0.359
        c_2 = -0.336

        BC = np.polynomial.polynomial.polyval(color_mag,[c_0,c_1,c_2])
              
        if np.logical_or(color_mag > 1.4,color_mag < -0.5).any():
            print "Warning : Out of range Color found"
        
        return BC
    
    if color == "g-r":  # valid range:(-0.2, 1.3). rms = 0.022
        
        c_0 = +0.053
        c_1 = -0.089
        c_2 = -0.736

        BC = np.polynomial.polynomial.polyval(color_mag,[c_0,c_1,c_2])
              
        if np.logical_or(color_mag > 1.3,color_mag < -0.2).any():
            print "Warning : Out of range Color found"
        
        return BC
    
    
    else:
        print "The first input is a string with the color used."
        return 0.    


def bol_correction_3(color,color_mag): # From Bersten-Hamuy
    '''From Bersten-Hamuy
    valid range B-V :(-0.2, 1.65). rms = 0.113'''

    if color == "B-V":  
        
        c_0 = -0.823
        c_1 = 5.027
        c_2 = -13.409
        c_3 = 20.133
        c_4 = -18.096
        c_5 = 9.084
        c_6 = -1.950
        
        BC = np.polynomial.polynomial.polyval(color_mag,[c_0,c_1,c_2,c_3,c_4,c_5,c_6])
              
        if np.logical_or(color_mag > 1.65,color_mag < -0.2).any():
            print "Warning : Out of range Color found"
        
        return BC
    else:
        print "The first input is a string with the color used."
        return 0.    
    
    
def Q_t(t,M_ni):
    '''
    t : Time in days (float or array)
    M_ni : Niquel mass in solar masses (float or array)
    Calculates the Luminosity from Ni+Co decay.
    Q(t) = M_ni*(6.45*np.exp(-t/t_ni)+1.45*np.exp(-t/t_co))*1e43  # [erg/s] . See Nadyozhin 1994
    returns 
    '''
    #day = 60*60*24 # day in seconds
    t_ni = 8.8        # niquel decay [days]
    t_co = 111.3      # Cobalt decay [days]
    t = np.asarray(t)
    Q = np.outer(np.asarray(M_ni),6.45*np.exp(-t/t_ni)+1.45*np.exp(-t/t_co))*1e43  # [erg/s] . See Nadyozhin 1994
    return Q 

def log_Q_t(t,M_ni):

    t_ni = 8.8        # niquel decay [days]
    t_co = 111.3      # Cobalt decay [days]
    Q = M_ni*(6.45*np.exp(-t/t_ni)+1.45*np.exp(-t/t_co))*1e43  # [erg/s] . See Nadyozhin 1994
    return np.log10(Q)

def M_bol_Q_t(t,M_ni):
    L_sun = 3.8270e33
    M_sun = 4.7554 
    return -2.5*(np.log10(Q_t(t,M_ni)/L_sun))+M_sun

def m_t(t,M_ni,T_0):
    
     return -2.5*(np.log10(Q_t(t,M_ni)/L_sun))+M_sun - 2.5*np.log10(1-np.exp(-(T_0/t)**2))
    

# LATE TIME MODEL    
    
def C_eta(eta):
    
    return (eta-3)**2. / (8*pi*(eta-1)*(eta-5))
 
def M_ejected(T_0,kappa=0.03,eta=10.,E=1e51): # T_0 = M * sqrt(C * kappa / E) 
    
    return T_0*np.sqrt(E/(C_eta(eta)*kappa))    
    
# CARDELLI LAW :                                                                                                                               


def F_a(x):

    if x<8 and x>5.9 :
        return -0.04473*(x-5.9)**2-0.009779*(x-5.9)**3
    else:
        return 0.0

def F_b(x):

    if x<8 and x>5.9 :

        return 0.2130*(x-5.9)**2+0.1207*(x-5.9)**3
    else:
        return 0.0
        
def a_x(x,verbose=False):

    if x > 0.3 and x <= 1.1:   # IR                                                                                                            

        return 0.574*x**1.61

    elif x > 1.1 and x <= 3.3:  # NIR/optical                                                                                                  

        y = x-1.82

        return 1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7

    elif x > 3.3 and x <=8 :   # UV                                                                                                            

        return 1.752 - 0.316*x-0.104/((x-4.67)**2+0.341) + F_a(x)

    else:

        if verbose:
            print "Out of range x"
        return 0.0



def b_x(x,verbose=False):

    if x > 0.3 and x <= 1.1:
        return -0.527*x**1.61

    elif x > 1.1 and x <= 3.3:

        y = x-1.82
        return 1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4 -0.62251*y**5+5.30260*y**6-2.09002*y**7

    elif x > 3.3 and x <=8 :

        return -3.090+1.825*x+1.206/((x-4.62)**2+0.263) +  F_b(x)
    else:

        if verbose:
            print "Out of range x"
        return 0.0


    
# ------------ -Black body --------------------------------------------- 

def B_lambda(lam,T): # cgs units 
    
    c = 2.99792458*1e10  # cm/s 
    k = 1.38064852*1e-16 
    h = 6.6260755e-27  # erg*s                             
    E = h*c/(lam*k*T)
    
    return np.pi * ((2. * h * c**2)/(lam**5))*(1/(np.exp(E)-1)) # erg/s/cm3


def B_nu(nu,T):
    
    c = 2.99792458*1e10
    k = 1.38064852*1e-16
    h = 6.6260755e-27
    E = h*nu/(k*T)
    
    return np.pi*((2.*h*nu**3)/c**2)*(1/(np.exp(E)-1))


def Black_body(lam,T,R,nu=False):
    
    #return (np.pi*R**2)*B_lambda(lam,T)
    if not nu:
        return (R**2)*B_lambda(lam,T)
    else:
        return (R**2)*B_nu(lam,T)


def int_BB(responses,T,R,nu=False,kernel=None,template=False,weight = 0.0):
    
    c = 2.99792458*1e10
    
    int_fluxes = []

    if not nu:
        
        for s_x,s_y in responses:     

            if s_x == None or s_y ==None:
                int_fluxes.append(0.0)
            else:
                BB_x = Black_body(s_x,T,R)
            
                if template:
                    kernel_x,kernel_y = kernel
                    kernel_X = interp1d(kernel_x,kernel_y)(s_x)
                    S_y = np.sqrt(BB_x * kernel_X) * s_y * s_x
                else:
                    S_y = BB_x * s_y * s_x
                
                I = integrate.trapz(S_y,s_x)
                int_fluxes.append(I)

        return np.asarray(int_fluxes)

    else :

        for s_x,s_y in responses:

            if s_x == None or s_y ==None:
                
                int_fluxes.append(0.0)
                
            else:
                s_nu = c / s_x
                BB_x = Black_body(s_nu,T,R,nu=nu)

                if template:
                    kernel_x,kernel_y = kernel
                    kernel_X = interp1d(kernel_x,kernel_y)(s_x)
                    #S_y = BB_x**(1-weight) * kernel_X**weight * s_y                                                                                                                                                                                                            
                    #S_y = BB_x * ((kernel_X/np.percentile(kernel_X,100)) * np.percentile(BB_x,100) / BB_x) **weight * s_y
                    S_y = BB_x * ((kernel_X/np.percentile(kernel_X,100))) **weight * s_y
                else:
                    S_y = BB_x * s_y

                I = integrate.trapz(S_y,np.log(s_nu))
                I_s = integrate.trapz(s_y,np.log(s_nu))
                int_fluxes.append(I/I_s)

        return np.asarray(int_fluxes)

    

# ---------  Late time parameters with free trapping ---------------




# Distance 

def power_law(x,a,b):

    return b*x**a

def power_law_2(x,a,b,c):

    return b*x**a+c

def E_B_V_poz(EW,R_V=3.1):
    
    E_B_V = 10**(2.16*EW-1.91)
    return E_B_V

def Photospheric_mag(t,t_0,filter='V'):

    if(filter=='V'):
        a_1 = -37.7  # (0.21)                                                                                                                                         \
                                                                                                                                                                       
        a_2 = 0.55  # (0.26)                                                                                                                                          \
                                                                                                                                                                       
        a_3 = -2.85  #(0.23)                                                                                                                                          \
                                                                                                                                                                       
    elif(filter=='I'):
            
        a_1 = -37.7 # (0.27)                                                                                                                                          \
                                                                                                                                                                       
        a_2 = 0.55  # (0.31)                                                                                                                                          \
                                                                                                                                                                       
        a_3 = -2.85  # (0.34)                                                                                                                                         \
                                                                                                                                                                       
    else:
        print "No calibration for %s filter"%(filter)
        return 0.

    M = a_1 + a_2*(t-t_0)/100. + a_3*np.log10((t-t_0)/100.)
    return M

def PMM_distance(m_corr,phot_vel,t,t_0,filter):  # m_corr : magnitude AKA corrected                                                                                    

    M = Photospheric_mag(t,t_0,filter)
    R = -5*np.log10(phot_vel)+37.76418
    mu = m_corr-R-M
    return mu

def inter_sample(x_0,x,y,yerr,n_sample = 1000,plot=False,confidence=True,n_c =95,honest_mean=True):
    
    inter = np.zeros(n_sample)
    sample_0 = np.interp(x_0,x,y)
    for i in np.arange(n_sample):
        
        y_sample = np.random.multivariate_normal(y,np.diag(yerr))
        inter_i = np.interp(x_0,x,y_sample)
        inter[i] = inter_i
        
    if plot:    
        pl.hist(inter,bins=10)
        pl.show()
    mean_inter = np.mean(inter)
    if confidence:
        std_inter = np.percentile(inter,q=[n_c,100-n_c])-np.mean(inter)
    else:
        std_inter = np.std(inter)
    if honest_mean:
        return sample_0,std_inter
    else:
        return mean_inter,std_inter
    
def fit_sample(x,y,yerr,model,n_par,n_sample = 1000,plot=False,confidence=True,n_c =95):
    
    inter = np.zeros((n_sample,n_par))
    
    for i in np.arange(n_sample):
        
        y_sample = np.random.multivariate_normal(y,np.diag(yerr))
        params,cov = curve_fit(f=linear_model,xdata=x,ydata=y_sample,sigma=yerr)
        inter[i] = params
        
    if plot:    
        for j in np.arange(n_par):
            pl.hist(inter[:,j],bins=int(n_sample/3.0))
            pl.show()
            pl.close()
    mean_inter = np.mean(inter,axis=0)
    if confidence:
        std_inter = np.percentile(inter,q=[n_c,100-n_c])-mean_inter
    else:
        std_inter = np.std(inter,axis=0)
    return mean_inter,std_inter

def model_interp(x0,x,y,yerr,model,n_sample = 1000,plot=False,maxfev = 1000):
    
    inter = np.zeros(n_sample)
    failed = 0
    if plot:
        figplot = pl.figure()
        axplot = figplot.add_subplot(111)
    for i in np.arange(n_sample):
        
        try:
            
            y_sample = np.random.multivariate_normal(y,np.diag(yerr))
            params,cov = curve_fit(f=linear_model,xdata=x,ydata=y_sample,sigma=yerr,maxfev = int(maxfev))
            inter[i] = linear_model(x0,*params)
            
        except:
            
            inter[i] = 0.0
            failed+=1
        
    if plot:
        
        axplot.hist(inter,bins=10)
        pl.show()
        
    print  "%s failed"%failed
    return np.sum(inter)/(len(inter)-failed),np.std(inter)

def model_interp_dif(x0,x1,x,y,yerr,model,n_sample = 1000,plot=False,maxfev = 1000):
    
    inter = np.zeros(n_sample)
    failed = 0
    if plot:
        figplot = pl.figure()
        axplot = figplot.add_subplot(111)
        figplot2 = pl.figure()
        axplot2 = figplot2.add_subplot(111)
    for i in np.arange(n_sample):
        
        try:
            
            y_sample = np.random.multivariate_normal(y,np.diag(yerr))
            params,cov = curve_fit(f=linear_model,xdata=x,ydata=y_sample,sigma=yerr,maxfev = int(maxfev))
            inter[i] = linear_model(x0,*params)-linear_model(x1,*params)
            if plot:
                baseline= np.arange(np.min([x0,x1]),np.max([x0,x1]),0.5)
                axplot2.plot(baseline,linear_model(baseline,*params))
            
        except:
            
            inter[i] = 0.0
            failed+=1
        
    if plot:
        
        axplot.hist(inter,bins=10)
        pl.show()
        
    print  "%s failed"%failed
    return np.sum(inter)/(len(inter)-failed),np.std(inter)

        
def fit_inter_sample(x0,x,y,yerr,model,n_sample = 1000,plot=False,maxfev = 1000):
    
    inter = np.zeros(n_sample)
    
    for i in np.arange(n_sample):
        
        y_sample = np.random.multivariate_normal(y,np.diag(yerr))
        params,cov = curve_fit(f=linear_model,xdata=x,ydata=y_sample,sigma=yerr,maxfev = int(maxfev))
        inter[i] = linear_model(x0,*params)
        
    if plot:    
        
        pl.hist(inter,bins=int(n_sample/3.0))
        pl.show()
        pl.close()
    mean_inter = np.mean(inter)
    std_inter = np.std(inter)
    return mean_inter,std_inter
    # Global variables 
    
def gauss_model(x,A,x0,fwhm):  # FWHM = 2sqrt(2ln(2))*sigma
    
    return A*np.exp(-4*np.log(2)*(x-x0)**2/fwhm)

def skew_gauss(x,A,x0,fwhm,alpha):
    
    def cum_gauss(x):
        
        return 0.5*(1+erf(x/np.sqrt(2)))
    
    sigma = fwhm / (2*np.sqrt(2*np.log(2.)))
    return 2.*gauss_model(x,A,x0,fwhm)*cum_gauss(alpha*(x-x0)/sigma)
    
    
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
