import io                           #Pipeline Module
import subprocess                   #Subprocess Module
import time                         #Timing Module
from PMOD_26 import *               #Parsing Module
  

# server set-up 
server = subprocess.Popen("./iso_server", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr = subprocess.STDOUT)

# Parameter Read-in

par_iso = list_file_grab('par_iso.don',[],False,True)

control_seq = par_iso[1]
control_iso = par_iso[4]

# Read-in Control values         
n_control = control_seq[0]       #
iphen_print = control_seq[1]     #
iso_calc = control_seq[2]        #
iso_print = control_seq[3]       # 
parab_print = control_seq[4]     #

# Read-in Iso values
n = control_iso[0]               #
mic = control_iso[1]             #
isnm = control_iso[2]            #
isym_emp = control_iso[3]        #
k0 = control_iso[4]              # 
rho0 = control_iso[5]            # 

# 
control_seq_str = str(n_control)+'\n'+str(iphen_print)+'\n'+str(iso_calc)+'\n'+str(iso_print)+'\n'+str(parab_print)+'\n'
control_iso_str = str(n)+'\n'+str(mic)+'\n'+str(isnm)+'\n'+str(isym_emp)+'\n'+str(k0)+'\n'+str(rho0)+'\n'

server.stdin.writelines([control_seq_str])
server.stdin.writelines([control_iso_str])  

#time_0 = time.time() 
#server.communicate("1\n0\n0\n0\n0\n0\n0\n") 
#time_1 = time.time() 







