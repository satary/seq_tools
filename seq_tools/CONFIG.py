import os,sys
EMAIL = "user@example.com" 

#Config assumes executable are in the same folder as python
#This usually happens if they are installed through conda.

# TEMP_DIR=os.path.expanduser('~/temp')
TEMP_DIR='/tmp'

MUSCLE_BIN_PATH=os.path.join(sys.exec_prefix,os.path.join('bin','muscle'))
PATH_TO_AL2CO=os.path.join(sys.exec_prefix,os.path.join('bin','al2co'))
BLOSSUM_PATH=os.path.expanduser('~/soft/al2co/BLOSUM62.txt')
