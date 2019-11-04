#!/usr/bin/env python3
import numpy as np
import sys
import os
from subprocess import DEVNULL, STDOUT, call, run, Popen

def which(program):
# This function checks for the existence of an executable
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return ''

def cleanup():
  for y in os.listdir("."):
    if os.path.isdir(y):
      if y!="correct_output" and y!="to_read":
        os.system("rm -r %s"%y)
    elif y!="test.inp" and y[-3:]!="dat":
      os.system("rm '%s'"%y)

print ("Welcome to the CONAN test script. Fingers crossed!")
for prog in ["gmx", "mencoder", "gnuplot"]:
  print ("Checking for %s..."%prog, end = '')
  if (which(prog) != ''):
    print("FOUND")
  else:
    print("MISSING")
    if (prog == "gmx"):
      print("Without gmx, CONAN can only work using the keyword REREAD. Testing it will not work now. Best of luck!")
      sys.exit()
    else:
      print("This is only used for visualization purposes. We can still test the program, but you are missing out!")

print ("Test 1: I/O, Result: ", end = '')
os.chdir("1_ReadWrite") # the first test, getting us the actual data.
run("../../conan.py test.inp".split(), stderr = STDOUT, stdout = DEVNULL, input = b"2\n")
i = run(["../compare_all.py"]).returncode
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
total_fails = 0
failures = []
print("Removing files from old tests...")
for test in tests:
  os.chdir(test)
  cleanup()
  os.chdir("..")
for test in tests:
  os.chdir(test)
  i += 1
  inpf = open("test.inp")
  header = inpf.readline()
  print ("Test %i: %s... "%(i, header[1:-1]), end = '')
  run("ln -s ../1_ReadWrite/dmf.xpm .".split())
  run("../../conan.py test.inp".split(), stderr = STDOUT, stdout = DEVNULL, input = b"2\n")
  e = run("python ../compare_all.py".split()).returncode
  if (e==0):
    print ("PASSED")
    os.system("rm dmf.xpm")
    cleanup()
  else:
    print ("FAILED")
    total_fails += 1
    failures.append(test)
  os.chdir("..")

if (total_fails == 0):
  print("All tests passed! Everything looks good.")
  os.chdir("1_ReadWrite/")
  cleanup()
else:
  print("%i test(s) failed. The output of the failed tests has been preserved."%len(failures))
  print("It is possible that the deviations are minor.")
  print("The failures are in the following folders:")
  for failure in failures:
    print(failure)
