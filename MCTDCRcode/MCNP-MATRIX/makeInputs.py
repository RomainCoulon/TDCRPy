import numpy as np
import linecache as lc

"""
Matrix basse energie
"""

#energy_i = np.linspace(1, 200, 200)        # emitted energy / keV
start_E = 0                              # start for the deposited energy / keV
end_E = 10000                               # end for the deposited energy / keV
delta_E= 10                             # the step the distribution / keV
#npas = int((end_E - start_E)/delta_E)
energy_i = np.linspace(2000, 10000, 801)        # emitted energy / keV
NPS = 1e6                                # number of source particlules for the MC calculation


for i in energy_i:
    with open("input3/input_"+str(int(i))+"keV.mcn", "w") as file:
        for j in range(1,44):
            a = lc.getline('template.txt',j)
            file.write(a)
        file.write("SI1 L "+str(round(i/1e3,3))+"             $ Energy in MeV\n") # line 44

        for l in range(45,49):
            a = lc.getline('template.txt',l) 
            file.write(a)
        file.write("E18:P,E "+str(start_E*1e-3)+" "+str(int(end_E/delta_E))+"i "+str((end_E+1)*1e-3)+"       $ binning of energy\n") # line 49
        file.write("NPS "+str(int(NPS))+"          $ number of source particles") # line 50
    
#f = open("output1/output_"+"keV.txt",'w')
    

with open("script3.bat", "w") as file2:
    for i in energy_i:
        file2.write("mcnp6 i=input3/input_"+str(int(i))+"keV.mcn r=run/run_"+str(int(i))+"keV.r o=output3/output_"+str(int(i))+"keV.o\n")