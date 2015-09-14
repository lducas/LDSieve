#!/bin/python

## Example usage: Benchmark.py .44 .44 > 44.bench
## stderr is quite verbose
## results of the benchmarks summarized in 44.bench

import subprocess
from math import *
import sys
print 

Targets = 100*[0]

Targets[50] = 3584095
Targets[51] = 3551524
Targets[52] = 3633605
Targets[53] = 3496843
Targets[54] = 3694084
Targets[55] = 3773021
Targets[56] = 3900625
Targets[57] = 3815991
Targets[58] = 4072324
Targets[59] = 3781187
Targets[60] = 3779136
Targets[61] = 4464769
Targets[62] = 4380649
Targets[63] = 4228565
Targets[64] = 4426816
Targets[65] = 4396757
Targets[66] = 4405628
Targets[67] = 4787344
Targets[68] = 4588164
Targets[69] = 4778537
Targets[70] = 4596736
Targets[71] = 4963938
Targets[72] = 4752400
Targets[73] = 4800481
Targets[74] = 5085025
Targets[75] = 5202961
Targets[76] = 5026564
Targets[77] = 5500000
Targets[78] = 5171076
Targets[79] = 5508409
Targets[80] = 5166529



a = float(sys.argv[1])
b = float(sys.argv[2])

spf = 4
d = 1./2.
m = 3
max_vec = 999999


def ttry(n,P,a,b):
	ss = 'cat ../challenges/dim'+str(n)+'sd0-LLL.txt | ../ldsieve ' 
	ss += ' ' +str(max_vec)
	ss += ' ' + str(Targets[n]) 
	ss += ' ' +str(P) 
	ss += ' ' +str(a) 
	ss += ' ' +str(b)
	print ss
	proc = subprocess.Popen(ss,shell = True, stdout=subprocess.PIPE)
	for line in proc.stdout:
		return [float(x) for x in line.split()]
	return ""


g = sqrt((a*a + b*b - 2. * a*b*d )/(1 - d*d) )
print "alpha = ",a
print "beta  = ",b
print "gamma = ",g

def C(n) :
	return sqrt(1-g*g)**n

print "[Cycle Cnt, time(sec), dim,  P  ,alph,beta, Nvecs , HashT Size, NbReduce ]"

for n in range(51,73,m):
	P = ceil(spf / (C(n)**(1./m)))
	print ttry(n,P,a,b)
	sys.stdout.flush()