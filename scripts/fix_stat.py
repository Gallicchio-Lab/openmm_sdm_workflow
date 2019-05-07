import os
import sys
import pickle
import re
import fileinput

#usage: python fix_stat.py < basename_stat.txt > basename_stat

stateid = {}
cycle = {}

status = {}

for line in fileinput.input():
    line.rsplit()
    items = line.split()
    try:
        k = int(items[0])
        stateid[k] = int(items[1])
        cycle[k] =  int(items[10])
    except:
        pass

nreplicas = len(cycle)
status = [{'stateid_current': stateid[k], 'running_status': 'W',
           'cycle_current': cycle[k]} for k in range(nreplicas)]

print(status)
pickle.dump(status,sys.stdout)
