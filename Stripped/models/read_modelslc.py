import numpy as np
import matplotlib.pyplot as pl

lcs = np.genfromtxt("lc.dat",dtype=[('file','S50')])
bands = ('U','B','V','R','I','J','H','K')
colormap = pl.cm.spectral
colors = [colormap(i) for i in np.linspace(0.1, 0.9,len(bands))]

for lc_file in lcs['file_']:
    
    print lc_file
    lc = np.genfromtxt(lc_file,dtype={'names':('file','t','Lbol','Luvoir','U','B','V','R','I','J','H','K'),\
                                          'formats':('S50','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8')})
    for band,color in zip(bands,colors):
        print band
        band_file = lc_file.split("_lc")[0]+'_%s.out'%band
        X=lc[['t','%s'%band]]
        np.savetxt(fname=band_file,X= np.append(zip(X['t'],X['%s'%band]),np.random.uniform(0.001,0.002,(X.shape[0],1)),axis=1),fmt=("%3.3f"))
        pl.plot(X['t'],X['%s'%band],label=band,marker='o',color=color)

    pl.legend()
    pl.gca().invert_yaxis()        
    pl.xlabel(r'$\mathrm{Time \ since \ explosion}$')
    pl.ylabel(r'$\mathrm{Magnitude \ (arbitrary)}$')
    pl.title(r'$%s$'%(lc_file.split("_lc")[0]))
    pl.savefig(lc_file.split("_lc")[0]+'_lcs.png')
    pl.show()
