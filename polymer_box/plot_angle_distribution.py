import mdtraj as md
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import numpy as np
from scipy.signal import lfilter,savgol_filter


def denoise(y):
    n = 15  # the larger n is, the smoother curve will be
    b = [1.0 / n] * n
    a = 1
    yy = lfilter(b, a, y)
    return yy

def denoise2(y):
    yy = savgol_filter(y, 31, 2)
    return yy


fins1 = [
'angle_1d-1nm-1.csv',
'angle_1d-1nm-2.csv',
'angle_1d-1nm-3.csv',
'angle_1d-1nm-4.csv',
'angle_1d-1nm-5.csv',
'angle_1d-1nm-6.csv',
'angle_1d-1nm-7.csv',
]

fins05 = [
'angle_1d-05nm-1.csv',
'angle_1d-05nm-2.csv',
'angle_1d-05nm-3.csv',
'angle_1d-05nm-4.csv',
'angle_1d-05nm-5.csv',
'angle_1d-05nm-6.csv',
'angle_1d-05nm-7.csv',
]
fins = fins05

data = []
for n,f in enumerate(fins):
    x1,y1,x2,y2 = np.loadtxt(f, delimiter=',',unpack=True)
    data.append( [f'poly-{n+1}',x1,y1,x2,y2] )

cmap = plt.get_cmap('jet')
colormap = cmap(np.linspace(0, 1, len(fins)))
color = [ colormap[c] for c in range(len(fins)) ] #'skyblue'
color = ['r','g','b','m','c','y','gray']

fig, axs = plt.subplots(1,2,figsize=(5,2.5))#,tight_layout=True,dpi=100)

base_line1,base_line2 = None, None
for n,d in enumerate(data):
    d[2] = d[2] / np.trapz( d[2], x=d[1] )
    d[4] = d[4] / np.trapz( d[4], x=d[3] )

    # Plane angles are 0-90
    mask = d[3]>90  # for planar orientation, angles are no larger than 90
    d[3] = np.abs( d[3]-mask*180 ) ## if angle>90, then use 180-angle so that all angles are 0<a<90
    angles = {}
    for x,y in zip(d[3],d[4]):
        if x in list( angles.keys() ):
            angles[x] += [y]
        else:
            angles[x] = [y]
    d[3] = list( angles.keys() )
    d[4] = [ np.mean(v) for v in angles.values() ]

    if base_line1 is None:
        base_line1 = np.copy(d[2])
    if base_line2 is None:
        base_line2 = np.copy(d[4])

    #axs[0].plot(d[1],d[2], color=color[n], lw=1, ls='-', marker='', ms=3, label=d[0], alpha=0.5)
    #axs[1].plot(d[3],d[4], color=color[n], lw=1, ls='-', marker='', ms=3, label=d[0], alpha=0.5)
    #axs[0].plot(d[1],d[2]-base_line1, color=color[n], lw=1, ls='-', marker='', ms=3, label=d[0])
    #axs[1].plot(d[3],d[4]-base_line2, color=color[n], lw=1, ls='-', marker='', ms=3, label=d[0])

    y1 = denoise2(d[2])
    y2 = denoise2(d[4])
    axs[0].plot(d[1],y1, color=color[n], lw=1, ls='-', marker='', ms=3, label=d[0])
    axs[1].plot(d[3],y2, color=color[n], lw=1, ls='-', marker='', ms=3, label=d[0])

    axs[0].set_xlabel('V$_A$-V$_A$ angle $(\degree)$',fontsize=10) ## input X name
    axs[1].set_xlabel('V$_B$-V$_B$ angle $(\degree)$',fontsize=10) ## input X name
    axs[0].set_ylabel('Distribution density',fontsize=10) ## input Y name
    axs[1].set_yticklabels([])
    for ii in [0,1]:
        axs[ii].tick_params(direction='in',labelsize='8' )
        axs[ii].set_ylim([0,0.012])
    axs[0].set_xlim([-7,187])
    axs[1].set_xlim([-7,90])
    axs[0].set_xticks( range(0,181,30) )
    axs[1].set_xticks( range(0,91,15) )


plt.subplots_adjust(wspace=0,bottom=0.15)#, hspace=0.0)

#axs[0].legend( loc='upper left', fontsize=10, frameon=False,)#  bbox_to_anchor=(0.99, 0.5))
plt.savefig('angle_distribution_1d.png', dpi=200)

exit()

print('Start 2D plot')
fins1 = [
'angle_2d-1nm-1.csv',
'angle_2d-1nm-2.csv',
'angle_2d-1nm-3.csv',
'angle_2d-1nm-4.csv',
'angle_2d-1nm-5.csv',
'angle_2d-1nm-6.csv',
'angle_2d-1nm-7.csv',
]
fins05 = [
'angle_2d-05nm-1.csv',
'angle_2d-05nm-2.csv',
'angle_2d-05nm-3.csv',
'angle_2d-05nm-4.csv',
'angle_2d-05nm-5.csv',
'angle_2d-05nm-6.csv',
'angle_2d-05nm-7.csv',
]
fins = fins05[:4]

fig, axs = plt.subplots(1,len(fins),figsize=(3*len(fins),3),tight_layout=True,dpi=100)

data = {}
base_line1 = None
for n,f in enumerate(fins):
    x1,y1 = np.loadtxt(f, delimiter=',',unpack=True)
    #data.append( [f'poly-{n+1}',x1,y1] )
    #axs[n].plot(xy[0],xy[1], color=color[n], lw=1, ls='', marker='o', ms=3, label=label[n], alpha=0.6)
    hist,xedge,yedge,_ = axs[n].hist2d(x1,y1, bins=180, density=True )
    if base_line1 is None:
        base_line1 = np.copy( hist )
    data[f] = [hist,xedge,yedge]
    #data[f] = [hist-base_line1,xedge,yedge]

plt.savefig('angle_distribution_2d.png', dpi=200)

import seaborn as sns

fig, axs = plt.subplots(1,len(data),figsize=(3*len(data),3),tight_layout=True,dpi=100)
for n,(k,v) in enumerate( data.items() ):
    sns.heatmap( v[0], ax=axs[n], cmap='viridis')
    #axs[n].imshow( v[0], cmap='viridis' )

plt.savefig('angle_distribution_2d_diff.png', dpi=200)

