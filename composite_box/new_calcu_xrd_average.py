import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from natsort import natsorted
from scipy import signal
from scipy import interpolate


def read_para(fin):
    with open(fin,'r') as ff1:
        lines = ff1.readlines()
        data = [ [] for n in range(len(lines[-1].split())) ]
        for line in lines[1:]:
            for n in range(len(data)):
                data[n].append( float(line.split()[n]) )
    return data


def write_para(x,y,fname):
  with open(fname, 'w') as ff2:
      #ff2.write("two_theta (degree)    intensity \n")
      for n in range(len(x)):
        ff2.write(f" {x[n]}  {y[n]} " + ' \n')


jdir = [
'job-box-poly1-5layeredGO/sep-gro-md-npt-03/',
#'job-box-poly3-5layeredGO/sep-gro-md-npt-03/',
'job-box-poly4-5layeredGO/sep-gro-md-npt-03/',
'job-box-poly1-5randomGO/sep-gro-md-npt-03/',
#'job-box-poly3-5randomGO/sep-gro-md-npt-03/',
'job-box-poly4-5randomGO/sep-gro-md-npt-03/',
]

xnew = np.linspace(2,78,num=1000) 

xrd = {} 
for j in jdir[:]:
    print( j )
    ## get two data first
    xrd_all, xrd_removed = [], []
    for x in os.listdir(j):
        if 'xrd-snapshot' in x:
            if 'removeBDOandSOL' in x:
                xrd_removed.append( x )
            else:
                xrd_all.append( x )
    xrd_all = natsorted(xrd_all)
    xrd_removed = natsorted(xrd_removed)

    #xrd_all = xrd_all[:10]
    #xrd_removed = xrd_removed[:10]

    print( len(xrd_all), len(xrd_removed) )

    # data 1
    data_xrd_all = []
    for x in xrd_all:
        fin = os.path.join(j,x)
        x1,y1 = np.loadtxt(fin, delimiter=',',unpack=True)
        interp_func = interpolate.interp1d(x1,y1)
        ynew = interp_func(xnew)
        data_xrd_all.append( ynew )

    #ynew = data_xrd_all[-1]
    ynew = np.mean( data_xrd_all, axis=0 )
    yerr = np.mean( np.abs( np.array(data_xrd_all)-ynew ), axis=0 )

    data_xrd_all = np.transpose( [xnew,ynew,yerr] )

    # data 2
    data_xrd_removed = []
    for x in xrd_removed:
        fin = os.path.join(j,x)
        x1,y1 = np.loadtxt(fin, delimiter=',',unpack=True)
        interp_func = interpolate.interp1d(x1,y1)
        ynew = interp_func(xnew)
        data_xrd_removed.append( ynew )

    #ynew = data_xrd_removed[-1]
    ynew = np.mean( data_xrd_removed, axis=0 )
    yerr = np.mean( np.abs( np.array(data_xrd_removed)-ynew ), axis=0 )

    data_xrd_removed = np.transpose( [xnew,ynew,yerr] )

    job = j.split('/')[0]
    np.savetxt(f"average_xrd_all_{job}.csv", data_xrd_all, delimiter=",")
    np.savetxt(f"average_xrd_removed_{job}.csv", data_xrd_removed, delimiter=",")
 
exit()

