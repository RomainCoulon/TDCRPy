import numpy as np
import linecache as lc

a = lc.getline('input/template.txt',5)
with open("input/input_test.txt",mode='w') as file:
    for i in range(1,44):
        a = lc.getline('input/template.txt',i)
        file.write(a)

    