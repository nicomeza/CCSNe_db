import os 

lcs = [l.split() for l in open('lc.standardsystem.sesn_allphot.dat','r')]
snlist = {}
scope = globals()
for lc in lcs:
    sn = lc[0]                                  
    if sn in snlist.keys():
        snlist[sn]['mag'].append(lc)
    else:
        snlist[sn] = {'mag':[]}
        snlist[sn]['mag'].append(lc) 

for key in snlist.keys():
    bands = {}
    lc_file = open("SN%s.dat"%key,'w')
    for mag in snlist[key]['mag']:
        band = mag[1]
        lc_file.write('\t'.join(mag[1:])+"\n")
        if bands.has_key(band):
            bands[band].append(mag[2:-1])
            scope['lc_%s'%band].write('\t'.join(mag[2:-1])+"\n")
        else:
            bands[band] = [mag[2:-1]]
            scope['lc_%s'%band] = open("SN%s_mags_%s.out"%(key,band),'w')
            scope['lc_%s'%band].write('\t'.join(mag[2:-1])+"\n")
    for band in bands:
        scope['lc_%s'%band].close()
        try:
            os.rename("SN%s_mags_%s.out"%(key,band),"../SN%s/Cfa_mags_%s.out"%(key,band)) 
            print key
        except:
            print "No dir SN%s/"%key
    lc_file.close()
    


## ----- Read invidual LCs ------------#
