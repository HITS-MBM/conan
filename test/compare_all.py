#!/usr/bin/env python3
import numpy as np
import sys
import os

def compare_halfprec(file_a, file_b):
    if (not os.path.exists(file_b)):
        return -1
    a = np.loadtxt(file_a)
    b = np.loadtxt(file_b)
    if (len(a) != len(b)):
      return -2 # wrong dimension
    compare_ab = np.isclose(a, b, atol = 1.1e-3, rtol = 1.1e-3)
    errors = len(np.where(compare_ab == False)[0])
    return errors

errors_list = []
for x in os.walk("correct_output/"):
    for y in x[2]:
        file_a = x[0] + "/" + y
        file_b = x[0][x[0].find("/")+1:]+"/" + y
        errors = compare_halfprec(file_a, file_b)
        if (errors != 0):
            errors_list.append(errors)
            if (errors==-1):
                print("This file doesn't exist at all: %s. This probably means that CONAN could not run at all."%file_b)
            if (errors==-2):
                print("This file has the wrong dimensions (!). This probably means that the CONAN test suite is out of date. Send an angry e-mail to the developers."%file_b)
            else:
                print("file_a: %s, file_b: %s, errors:%i"%(file_a, file_b, errors))
sys.exit(len(errors_list))
