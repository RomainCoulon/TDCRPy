import numpy as np

"""
Matrix basse energie
"""
energy_i = np.linspace(0, 20, 20)        # emitted energy / keV
start_E = 0                              # start for the deposited energy / keV
end_E = 20                               # end for the deposited energy / keV
delta_E= 0.2                             # the step the distribution / keV
NPS = 1e5                                # number of source particlules for the MC calculation


for i in energy_i:
    with open("input/input_"+str(int(i))+"keV.txt", "w") as file:
        file.write("MatrixGamma\n")
        file.write("SI1 L "+str(round(i/1e3,3))+"             $ Energy in MeV\n") # line 44
        file.write("E18:P,E "+str(start_E)+" "+str(int(end_E/delta_E))+"i "+str(end_E)+"       $ binning of energy\n") # line 49
        file.write("NPS "+str(int(NPS))+"          $ number of source particles") # line 49


with open("script.bat", "w") as file2:
    for i in energy_i:
        file2.write("mcnp6 i=input/input_"+str(int(i))+"keV.mcn r=run/run_"+str(int(i))+"keV.r o=output/output_"+str(int(i))+"keV.o\n")