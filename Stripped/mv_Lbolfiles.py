import os 
import subprocess

sn_labels = ['sn','type','host','hostredshift','hostlumdist','sn_ebv','sn_z','t_0','t_discov','m_discov','t_non_det','m_non_det']
sn_formats = ['S15','S10','S20','f8','f8','f8','f8','f8','f8','f8','f8','f8']

SN_DATA = np.genfromtxt(IN_FILE,dtype={'names':sn_labels,'formats':sn_formats})

for SN,z_SN,E_B_V,t_0,d_L,t_discov,t_non in SN_DATA[['sn','sn_z','sn_ebv','t_0','hostlumdist','t_discov','t_non_det']]:

    print "####### %s ########## \n"%SN
    try:

        os.chdir("%s/"%SN)
        last_file = subprocess.check_output("ls -ltr | grep dat | tail -1 | awk '{print $NF}'", shell=True).split()[0]
        print last_file
        s = last_file.split("_")
        new_file = s[0]+"_Lbol_BVRIYJH.dat"
        os.system("cp %s %s"%(last_file,new_file))
        os.chdir(current_dir)
    except:
        os.chdir(current_dir)
                    
