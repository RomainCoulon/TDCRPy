import numpy as np
import linecache as lc

"""
Matrix basse energie
"""
'''
#energy_i = np.linspace(1, 200, 200)        # emitted energy / keV
start_E = 0                              # start for the deposited energy / keV
end_E = 10000                               # end for the deposited energy / keV
delta_E= 10                             # the step the distribution / keV
#npas = int((end_E - start_E)/delta_E)
energy_i = np.linspace(2000, 10000, 801)        # emitted energy / keV
NPS = 1e5                                # number of source particlules for the MC calculation


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
    
#f = open("output/output_"+"keV.txt",'w')
    

with open("script3.bat", "w") as file2:
    for i in energy_i:
        file2.write("mcnp6 i=input3/input_"+str(int(i))+"keV.mcn r=run/run_"+str(int(i))+"keV.r o=output3/output_"+str(int(i))+"keV.o\n")
'''

def make_input(niveau,NPS,par):
    ## par (string) : 'p' (photon) ou 'b' (beta-) ou 'bp' (beta+)
    start_E = 0

    if niveau == 0:
        end_energy = 200
        start_energy = 1
        taille_x = 200      # 1k-200k où delta=1k
        input_n = '1'
        delta_E = 0.2

    elif niveau == 1:
        end_energy = 2000
        start_energy = 200
        taille_x = 901     # 200k-2000k où delta=1k
        input_n = '2'
        delta_E = 2

    else:
        end_energy = 1e4
        start_energy = 2e3
        taille_x = 801      #2M-10M où delta = 0.1M
        input_n = '3'
        delta_E = 10

    energy_i = np.linspace(start_energy,end_energy,taille_x)

    if par == 'p':
        doc = 'template.txt'
    elif par == 'b':
        doc = 'template_b.txt'
    elif par == 'bp':
        doc = 'template_bp.txt'
    
    name_bat = "script" + input_n + par +".bat"

    for i in energy_i:
        name_doc = "input" + input_n + par + "/input_" + str(int(i)) + "keV.mcn"

        with open(name_doc, "w") as file1:
            for j in range(1,44):
                a = lc.getline(doc,j)
                file1.write(a)
            file1.write("SI1 L "+str(round(i/1e3,3))+"             $ Energy in MeV\n") # line 44

            for l in range(45,49):
                a = lc.getline(doc,l) 
                file1.write(a)
            file1.write("E18:P,E "+str(start_E*1e-3)+" "+str(int(end_energy/delta_E))+"i "+str((end_energy+1)*1e-3)+"       $ binning of energy\n") # line 49
            file1.write("NPS "+str(int(NPS))+"          $ number of source particles") # line 50

    with open(name_bat, "w") as file2:
        for i in energy_i:
            contenu = "mcnp6 i=input" + input_n + par + "/input_" + str(int(i)) + "keV.mcn r=run/run_"+str(int(i))+"keV.r o=output" + input_n + par + "/output_"+str(int(i))+"keV.o\n"
            file2.write(contenu)
    return 1

write = make_input(2,1e6,'bp')   
 