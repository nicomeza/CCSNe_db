import json as js
import numpy as np
import matplotlib.pyplot as pl
import os 

snlist = [f.split()[0][:-1] for f in open('sn_dirs.dat','r').readlines()]
homedir = "/home/dust_speck/SN/CCSNe_db/Stripped/"
#sn_data = open('sn_data.dat','w')
params = ['sn','host','hostredshift','hostlumdist','host_ebv']
#sn_data.write("%s \n"%"\t".join(params))
uvot_bands=['U','B','V','M2','W1','W2','UVW1','UVW2','UVM2','u','b','v']
no_bands = ['F475W','F625W','F775W','C','G']

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

        try:
            
            host_str = ','.join([x['value'] for x in json_sn['%s'%sn]['host']])
            host_z = ','.join([x['value'] for x in json_sn['%s'%sn]['hostredshift']])
            try:
                host_dL = ','.join(["%s [%s]"%(x['value'],x['u_value']) for x in json_sn['%s'%sn]['hostlumdist']])
            except:
                host_dL = ','.join(["%s [%s]"%(x['value'],x['u_value']) for x in json_sn['%s'%sn]['lumdist']])
            host_ebv = ','.join([x['value'] for x in json_sn['%s'%sn]['ebv']])
            print sn,host_str,host_z,host_dL,host_ebv
        
        except:
            print "no host information ?"
            
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
            system = phot.get("System")
            if band in uvot_bands:
                if teles == "Swift":
                    
                    if band == "UVW1":
                        band = "W1"
                    elif band == "UVW2":
                        band = "W2"
                    elif band == "UVM2":
                        band = "M2"
                    
                    band = "%s_uvot"%band
                    
            if band != None and time!= None and not (band in no_bands):
            
                if band in bands:
                    if phot['u_time'] != "MJD":
                        print ("NO MJD Here")
                    try:
                        bands[band].append((phot['time'],phot['magnitude'],phot['e_magnitude']))
                    except:
                        print "No error!"
                        bands[band].append((phot['time'],phot['magnitude'],"999"))
                else:
                    
                    bands["%s"%band] = []
                    print "adding %s"%band
                    
                    try:
                        bands[band].append((phot['time'],phot['magnitude'],phot['e_magnitude']))
                    except:
                        print ("No error!")
                        bands[band].append((phot['time'],phot['magnitude'],"999"))
            else:
            
                print "Warning : No band"
                print phot
                
        
        os.system("rm mags_*.out")
        
        fig = pl.figure(figsize=(10,6))
        ax = fig.add_subplot(111)
        colormap = pl.cm.spectral
        fig.gca().set_color_cycle([colormap(i) for i in np.linspace(0.1, 0.9,len(bands.keys()))])


        for key in iter(bands):
        
            print key
            try:
                t,mag,magerr = np.asfarray(bands.get(key))[:,0],np.asfarray(bands.get(key))[:,1],np.asfarray(bands.get(key))[:,2]
                with_err = np.where(magerr!=999)[0]
                no_err = np.where(magerr==999)[0]
                if len(with_err)>0:
                    ax.errorbar(t[with_err],mag[with_err],yerr=magerr[with_err],marker='o',label=key)
                if len(no_err)>0:
                    ax.errorbar(t[no_err],mag[no_err],marker='o',label="no err %s"%key,alpha=0.3)
                np.savetxt('mags_%s.out'%key, zip(t,mag,magerr),fmt="%5.3f")
                print 'Saved mags_%s.out'%key
                #print mag
            except:
                print "key : %s"%key
                t,mag = np.asfarray(bands.get(key))[:,0],np.asfarray(bands.get(key))[:,1]    
                ax.errorbar(t,mag,marker='o',label=key)
                np.savetxt('mags_%s.out'%key, zip(t,mag),fmt='%5.3f')
                print 'Saved mags_%s.out'%key
                
        ax.set_ylabel('mag')
        ax.set_xlabel('MJD')
        ax.set_title("%s"%sn)
        fig.gca().invert_yaxis()
        ax.legend(loc='best',ncol=2,prop={'size':8})
        os.system("rm %s_phot.png"%sn)
        pl.savefig("%s_phot.png"%sn)
        pl.close()
        os.chdir(homedir)

except:
    
    print "FATAL ERROR"
    os.chdir(homedir)
    pl.close()
#    sn_data.close()

#sn_data.close()
