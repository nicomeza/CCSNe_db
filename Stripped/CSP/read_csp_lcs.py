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