# Get a test case from Fortran GMET repo
# https://github.com/NCAR/GMET/tree/master/test_cases

import os

targetpath = '.'
os.chdir(targetpath)

os.system('git clone git@github.com:NCAR/GMET.git')
os.system('tar -xf  GMET/test_cases/cali2017.tgz')
os.system('rm -rf GMET')

