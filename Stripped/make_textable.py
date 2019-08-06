import numpy as np

IN_FILE = "sn_data.dat"

sn_labels = ['sn','type','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_head = ['SN','Type','Host','Host redshift','host $d_L$','$E(B-V)$','$t_0$','t-discovery','mag-discovery','t-non detection','m-non detection']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8']

description = "Sample of SNe"

SN_DATA = np.genfromtxt(IN_FILE,dtype={'names':sn_labels,'formats':sn_formats})

sn_tex = open('sne_table.tex','w')

columns = len(sn_labels)*'c'
colhead = '&'.join(['{%s}'%s for s in sn_head])

sn_tex.write("""\\begin{table*}\n\\centering\n\\scriptsize\n\\caption{%s\\label{tab:table}}\n\\begin{tabular}{%s}\n\\hline\n%s\n"""\
                 %(description,columns,colhead))

for sn_data in SN_DATA[['sn','type','host','hostredshift','hostlumdist','sn_ebv','t_0','t_discov','m_discov','t_non_det','m_non_det']]:
    
    data_s = '&'.join(['%s'%s for s in sn_data])
    sn_tex.write(data_s+"\\\ \n")



sn_tex.write("""\end{tabular}\n
\\tablenotetext{a}{Heliocentric Redshift}\n
\\tablecomments{Comments}\n\end{table}\n""")
               
sn_tex.close()
