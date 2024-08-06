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


path_pbi = '/qfs/projects/sepcon/difan/sepcon-pbi/work-gmx/jobs-poly/'

fins05 = [
path_pbi+'angle_1d-05nm-1.csv',
'angle_1d-05nm-poly1-5randomGO.csv',
'angle_1d-05nm-poly1-5layeredGO.csv',

#path_pbi+'angle_1d-05nm-4.csv',
#'angle_1d-05nm-poly4-5layeredGO.csv',
#'angle_1d-05nm-poly4-5randomGO.csv',
]
fins = fins05

data = []
for n,f in enumerate(fins):
    x1,y1,x2,y2 = np.loadtxt(f, delimiter=',',unpack=True)
    data.append( [f'poly-{n+1}',x1,y1,x2,y2] )
data[2][2] = data[2][2]*0.9 + data[1][2]*0.1

cmap = plt.get_cmap('jet')
colormap = cmap(np.linspace(0, 1, len(fins)))
color = [ colormap[c] for c in range(len(fins)) ] #'skyblue'
color = ['r','g','b','m','c','y','gray']
color = ['r','seagreen','skyblue']

fig, axs = plt.subplots(1,2,figsize=(3,2))#,tight_layout=True,dpi=100)

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

    axs[0].set_xlabel('                 Angles $(\degree)$',fontsize=10) ## input X name
    #axs[1].set_xlabel('V$_B$-V$_B$ \nangle $(\degree)$',fontsize=10) ## input X name
    axs[0].set_ylabel('Angle distribution',fontsize=10) ## input Y name
    axs[1].set_yticklabels([])
    for ii in [0,1]:
        axs[ii].tick_params(direction='in',labelsize='8' )
        axs[ii].set_ylim([0,0.012])
    axs[0].set_xlim([-5,200])
    axs[1].set_xlim([-10,90])
    axs[0].set_xticks( [0,60,120,180] )
    axs[1].set_xticks( [0,30,60,90] )


plt.subplots_adjust(wspace=0,bottom=0.2, left=0.2)#, hspace=0.0)

#axs[0].legend( loc='upper left', fontsize=10, frameon=False,)#  bbox_to_anchor=(0.99, 0.5))
plt.savefig('angle_distribution_1d.png', dpi=800)

