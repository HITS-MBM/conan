#!/usr/bin/env python
import numpy as np
import sys
import os

print ("Welcome to the CONAN test script. Fingers crossed!")
print ("Test 1: I/O, Result: ", end = '')
os.chdir("1_ReadWrite") # the first test, getting us the actual data.
os.system("echo 2 | ../../conan.py test.inp > /dev/null >& /dev/null")
i = os.system("python ../compare_all.py ")
if (i==0):
  print ("PASSED")
else:
  print ("FAILED")

i = 1
tests = []
os.chdir("..")
for x in os.listdir("."):
  if os.path.isdir(x) and x!="1_ReadWrite":
    tests.append(x)
tests.sort()
for test in tests:
  os.chdir(test)
  i += 1
  inpf = open("test.inp")
  header = inpf.readline()
  print ("Test %i: %s... "%(i, header[1:-1]), end = '')
  os.system("ln -s ../1_ReadWrite/dmf.xpm .")
  os.system("../../conan.py test.inp > /dev/null >& /dev/null")
  e = os.system("python ../compare_all.py")
  if (e==0):
    print ("PASSED")
  else:
    print ("FAILED")
  os.system("rm dmf.xpm")
  for y in os.listdir("."):
    if os.path.isdir(y) and y!="correct_output" and y!="to_read":
      os.system("rm -r %s"%y)
  os.chdir("..")
