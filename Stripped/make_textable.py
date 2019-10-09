import numpy as np

IN_FILE = "sn_test_total.dat"

sn_labels = ['sn','type','host','hostredshift','hostlumdist','MW_ebv','host_ebv','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_head = ['SN','Type','Host','Host redshift','host $d_L$','$MW E(B-V)$','$host E(B-V)$','$t_0$','t-discovery','mag-discovery','t-non detection','m-non detection']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8']

description = "Sample of SNe"

SN_DATA = np.genfromtxt(IN_FILE,dtype={'names':sn_labels,'formats':sn_formats})

sn_tex = open('sne_table.tex','w')

columns = len(sn_labels)*'c'
colhead = '&'.join(['{%s}'%s for s in sn_head])

sn_tex.write("""\\begin{table*}\n\\centering\n\\scriptsize\n\\caption{%s\\label{tab:table}}\n\\begin{tabular}{%s}\n\\hline\n%s\n"""\
                 %(description,columns,colhead))

for sn_data in SN_DATA[['sn','type','host','hostredshift','hostlumdist','MW_ebv','host_ebv','t_0','t_discov','t_non_det']]:
    
    sn,sn_type,host,host_z,host_dL,mw_ebv,host_ebv,t_0,t_disc,t_nondet = sn_data

    if t_0<100.:
        t_0 = 0.5*(t_disc+t_nondet)
    host = host[1:-1]
    data_s = ' & '.join(['%s'%s for s in [sn,sn_type,host,host_z]])+' & '
    data_s += ' & '.join(['%2.2f'%s for s in [host_dL,mw_ebv,host_ebv,t_0]])
    sn_tex.write(data_s+"\\\ \n")



sn_tex.write("""\end{tabular}\n
\\tablenotetext{a}{Heliocentric Redshift}\n
\\tablecomments{Comments}\n\end{table}\n""")
               
sn_tex.close()
