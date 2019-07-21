lcs = [l.split() for l in open('table8.dat','r')]
snlist = {}

for lc in lcs:
    sn = lc[0]                                  
    if sn in snlist.keys():
        snlist[sn]['mag'].append(lc)
    else:
        snlist[sn] = {'mag':[]}
        snlist[sn]['mag'].append(lc) 

for key in snlist.keys():
    lc_file = open("SN%s_NIR_duPont.dat"%key,'w')
    for mag in snlist[key]['mag']:
        lc_file.write('\t'.join(mag[1:])+"\n")
    lc_file.close()
    

#  -------------- NIR du pont -------------------------------

lcs = [l.split()[0] for l in open("NIR_duPont.dat")]

bands = ["Y","J","H"]
for file in lcs:
    
    sn = file.split("_")[0]
    print sn
    mags = np.genfromtxt("%s_NIR_duPont.dat"%sn,dtype=[("JD","f8"),("Y","f8"),("Yerr","f8"),("J","f8"),("Jerr","f8"),("H","f8"),("Herr","f8")])
    mags['JD']-= 2400000.5
    if not os.path.exists("../../%s"%sn):
        os.makedirs("../../%s"%sn)
    for band in bands:
        try:
            data = mags[['JD',"%s"%band,"%serr"%band]]
            where_dat = np.where(~np.isnan(mags["%s"%band]))
            np.savetxt("%s_mags_%s.out"%(sn,band),data[where_dat])
            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_WIRC.out"%(sn,band)) 
        except:
            try :
                jd,mag,err = data['JD'],data['%s'%band],data["%serr"%band]
                if mag!=np.nan:
                    np.savetxt("%s_mags_%s.out"%(sn,band),np.array([[jd],[mag],[err]]).T)
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_WIRC.out"%(sn,band)) 
            except:    
                print "no data in band %s"%band

#  -------------- OPT du pont -------------------------------

lcs = [l.split()[0] for l in open("OPT_duPont.dat")]

bands = ["u","g","r","i","B","V"]
for file in lcs:
    
    sn = file.split("_")[0]
    print sn
    mags = np.genfromtxt("%s_OP_duPont.dat"%sn,\
                             dtype=[("JD","f8"),("u","f8"),("uerr","f8"),("g","f8"),("gerr","f8"),("r","f8"),("rerr","f8"),("i","f8"),("ierr","f8"),("B","f8"),("Berr","f8"),("V","f8"),("Verr","f8")])
    mags['JD']-= 2400000.5
    if not os.path.exists("../../%s"%sn):
        os.makedirs("../../%s"%sn)
    for band in bands:
        try:
            data = mags[['JD',"%s"%band,"%serr"%band]]
            np.savetxt("%s_mags_%s.out"%(sn,band),data)
            if band != "V":
                os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_CSPI.out"%(sn,band)) 
            else:
                if np.min(data['JD']<=53748.0):
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_3014.out"%(sn)) 
                elif np.min(data['JD']<=53759.0):
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_3009.out"%(sn)) 
                else:
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_9844.out"%(sn)) 
        except:
            try :
                jd,mag,err = data['JD'],data['%s'%band],data["%serr"%band]
                if mag!=np.nan:
                    
                    np.savetxt("%s_mags_%s.out"%(sn,band),np.array([[jd],[mag],[err]]).T)
                    
                    if band != "V":
                        os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_CSPI.out"%(sn,band)) 
                    else:
                        if np.min(data['JD']<=53748.0):
                            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_3014.out"%(sn)) 
                        elif np.min(data['JD']<=53759.0):
                            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_3009.out"%(sn)) 
                        else:
                            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_9844.out"%(sn)) 
                            
                            
            except:    
                print "no data in band %s"%band


#  -------------- NIR Swope -------------------------------

lcs = [l.split()[0] for l in open("NIR_Swope.dat")]

bands = ["Y","J","H"]
for file in lcs:
    
    sn = file.split("_")[0]
    print sn
    mags = np.genfromtxt("%s_NIR_Swope.dat"%sn,dtype=[("JD","f8"),("Y","f8"),("Yerr","f8"),("J","f8"),("Jerr","f8"),("H","f8"),("Herr","f8")])
    mags['JD']-= 2400000.5
    if not os.path.exists("../../%s"%sn):
        os.makedirs("../../%s"%sn)
    for band in bands:
        try:
            data = mags[['JD',"%s"%band,"%serr"%band]]
            np.savetxt("%s_mags_%s.out"%(sn,band),data)
            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/CSP_Swope_mags_%s.out"%(sn,band)) 
            if band == "Y" or band =="H":
                os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_RC.out"%(sn,band)) 
            else:
                if np.min(jd)<=56007.0:
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_RC1.out"%(sn,band)) 
                else:
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_RC2.out"%(sn,band)) 
        except:
            try :
                jd,mag,err = data['JD'],data['%s'%band],data["%serr"%band]
                if mag!=np.nan:
                    np.savetxt("%s_mags_%s.out"%(sn,band),np.array([[jd],[mag],[err]]).T)
                    if band == "Y" or band =="H":
                        os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_RC.out"%(sn,band)) 
                    else:
                        if np.min(jd)<=56007.0:
                            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_RC1.out"%(sn,band)) 
                        else:
                            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_RC2.out"%(sn,band)) 
            except:  
                print "no data in band %s"%band

#  -------------- OPT Swope  -------------------------------

lcs = [l.split()[0] for l in open("OPT_Swope.dat")]

bands = ["u","g","r","i","B","V"]
for file in lcs:
    
    sn = file.split("_")[0]
    print sn
    mags = np.genfromtxt("%s_OP_Swope.dat"%sn,\
                             dtype=[("JD","f8"),("u","f8"),("uerr","f8"),("g","f8"),("gerr","f8"),("r","f8"),("rerr","f8"),("i","f8"),("ierr","f8"),("B","f8"),("Berr","f8"),("V","f8"),("Verr","f8")])
    mags['JD']-= 2400000.5
    if not os.path.exists("../../%s"%sn):
        os.makedirs("../../%s"%sn)
    for band in bands:
        try:
            data = mags[['JD',"%s"%band,"%serr"%band]]
            np.savetxt("%s_mags_%s.out"%(sn,band),data)
            if band != "V":
                os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_CSPI.out"%(sn,band)) 
            else:
                if np.min(data['JD']<=53748.0):
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_3014.out"%(sn)) 
                elif np.min(data['JD']<=53759.0):
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_3009.out"%(sn)) 
                else:
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_9844.out"%(sn)) 

        except:
            try :
                jd,mag,err = data['JD'],data['%s'%band],data["%serr"%band]
                if mag!=np.nan:
                    np.savetxt("%s_mags_%s.out"%(sn,band),np.array([[jd],[mag],[err]]).T)
                    if band != "V":
                        os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_%s_CSPI.out"%(sn,band)) 
                    else:
                        if np.min(data['JD']<=53748.0):
                            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_3014.out"%(sn)) 
                        elif np.min(data['JD']<=53759.0):
                            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_3009.out"%(sn)) 
                        else:
                            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/mags_V_9844.out"%(sn)) 

            except:    
                print "no data in band %s"%band
