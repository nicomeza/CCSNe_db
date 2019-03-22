import json as js
import numpy as np
import matplotlib.pyplot as pl
import os 

snlist = [f.split()[0][:-1] for f in open('sndirs.dat','r').readlines()]
homedir = "/home/dustspeck/SN/CCSNe_db/Stripped/"
sn_data = open('sn_data.dat','w')
params = ['sn','host','hostredshift','hostlumdist','host_ebv']
sn_data.write("%s \n"%"\t".join(params))
uvot_bands=["U",'B','V']

try :
    for sn in snlist:
    
        print sn
        os.chdir("%s/"%sn)
        json_sn = js.loads(open('%s.json'%sn).read())
        
        try:

            phot_list = json_sn['%s'%sn]['photometry']
        except:
            print "Wrong SN name or no photometry ?"
            phot_list = json_sn['%s'%sn[2:]]['photometry']
            sn = sn[2:]

        host_str = ','.join([x['value'] for x in json_sn['%s'%sn]['host']])
        host_z = ','.join([x['value'] for x in json_sn['%s'%sn]['hostredshift']])
        try:
            host_dL = ','.join(["%s [%s]"%(x['value'],x['u_value']) for x in json_sn['%s'%sn]['hostlumdist']])
        except:
            host_dL = ','.join(["%s [%s]"%(x['value'],x['u_value']) for x in json_sn['%s'%sn]['lumdist']])
            
        host_ebv = ','.join([x['value'] for x in json_sn['%s'%sn]['ebv']])
        print sn,host_str,host_z,host_dL,host_ebv
        #sn_data.write("%s \n"%("\t".join((sn,host_str,host_z,host_dL,host_ebv))))

        try:

            phot_list = json_sn['%s'%sn]['photometry']
        except:
            print "Wrong SN name or no photometry ?"
            phot_list = json_sn['%s'%sn[2:]]['photometry']
        bands = {}

        for phot in phot_list:
        
        
            band = phot.get('band')
            time = phot.get('time')
            teles = phot.get('telescope')
            
            if band in uvot_bands:
                if teles == "Swift":
                    band = "%s_uvot"%band
                    
            if band != None and time!= None:
            
                if band in bands:
                    if phot['u_time'] != "MJD":
                        print "NO MJD Here"
                    try:
                        bands[band].append((phot['time'],phot['magnitude'],phot['e_magnitude']))
                    except:
                        bands[band].append((phot['time'],phot['magnitude'],"999"))
                else:
                    
                    bands["%s"%band] = []
                    print "adding %s"%band
                    
                    try:
                        bands[band].append((phot['time'],phot['magnitude'],phot['e_magnitude']))
                    except:
                        bands[band].append((phot['time'],phot['magnitude'],"999"))
            else:
            
                print "Warning : No band"
                print phot
                
        
        for key in iter(bands):
    
            print key
            try:
                t,mag,magerr = np.asfarray(bands.get(key))[:,0],np.asfarray(bands.get(key))[:,1],np.asfarray(bands.get(key))[:,2]
                with_err = np.where[magerr!=999][0]
                pl.errorbar(t[with_err],mag[with_err],yerr=magerr[with_err],marker='o',label=key)
                np.savetxt('mags_%s.out'%key, zip(t,mag,magererr),fmt="%5.3f") 
            except:
                t,mag = np.asfarray(bands.get(key))[:,0],np.asfarray(bands.get(key))[:,1]    
                pl.errorbar(t,mag,marker='o',label=key)
                np.savetxt('mags_%s.out'%key, zip(t,mag),fmt='%5.3f')
                
        
        pl.ylabel('mag')
        pl.xlabel('MJD')
        pl.title("%s"%sn)
        pl.gca().invert_yaxis()
        pl.legend()
        os.system("rm %s_phot.png"%sn)
        pl.savefig("%s_phot.png"%sn)
        #pl.show()
        os.chdir(homedir)

except:
    
    print "FATAL ERROR"
    os.chdir(homedir)
    pl.close()
    sn_data.close()

sn_data.close()