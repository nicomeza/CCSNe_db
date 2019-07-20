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
            np.savetxt("%s_mags_%s.out"%(sn,band),data)
            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/CSP_duPont_mags_%s.out"%(sn,band)) 
        except:
            try :
                jd,mag,err = data['JD'],data['%s'%band],data["%serr"%band]
                if mag!=np.nan:
                    np.savetxt("%s_mags_%s.out"%(sn,band),np.array([[jd],[mag],[err]]).T)
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/CSP_duPont_mags_%s.out"%(sn,band)) 
            except:    
                print "no data in band %s"%band

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
            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/CSP_duPont_mags_%s.out"%(sn,band)) 
        except:
            try :
                jd,mag,err = data['JD'],data['%s'%band],data["%serr"%band]
                if mag!=np.nan:
                    np.savetxt("%s_mags_%s.out"%(sn,band),np.array([[jd],[mag],[err]]).T)
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/CSP_duPont_mags_%s.out"%(sn,band)) 
            except:    
                print "no data in band %s"%band



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
        except:
            try :
                jd,mag,err = data['JD'],data['%s'%band],data["%serr"%band]
                if mag!=np.nan:
                    np.savetxt("%s_mags_%s.out"%(sn,band),np.array([[jd],[mag],[err]]).T)
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/CSP_Swope_mags_%s.out"%(sn,band)) 
            except:    
                print "no data in band %s"%band

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
            os.rename("%s_mags_%s.out"%(sn,band),"../../%s/CSP_Swope_mags_%s.out"%(sn,band)) 
        except:
            try :
                jd,mag,err = data['JD'],data['%s'%band],data["%serr"%band]
                if mag!=np.nan:
                    np.savetxt("%s_mags_%s.out"%(sn,band),np.array([[jd],[mag],[err]]).T)
                    os.rename("%s_mags_%s.out"%(sn,band),"../../%s/CSP_Swope_mags_%s.out"%(sn,band)) 
            except:    
                print "no data in band %s"%band
