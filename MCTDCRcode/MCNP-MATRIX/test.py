import numpy as np
import linecache as lc

energy_i = np.linspace(0, 20, 20)        # emitted energy / keV
start_E = 0                              # start for the deposited energy / keV
end_E = 20                               # end for the deposited energy / keV
delta_E= 0.2                             # the step the distribution / keV
NPS = 1e5 


with open("input/input_test.txt",mode='w') as file:
    for i in range(1,44):
        a = lc.getline('input/template.txt',i)
        file.write(a)
    file.write("SI1 L "+str(15)+"             $ Energy in MeV\n") # line 44

    for i in range(45,49):
        a = lc.getline('input/template.txt',i) 
        file.write(a)
    file.write("E18:P,E "+str(15)+" "+str(int(15))+"i "+str(20)+"       $ binning of energy\n") # line 49
    file.write("NPS "+str(int(15))+"          $ number of source particles") # line 50


    