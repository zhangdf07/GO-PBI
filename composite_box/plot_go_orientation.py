import mdtraj as md
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import numpy as np
from scipy.signal import lfilter,savgol_filter


def denoise2(y):
    yy = savgol_filter(y, 11, 2)
    return yy


fins = [
'go_orientation-poly1-5randomGO.csv',
'go_orientation-poly1-5layeredGO.csv',
]

data = []
for n,f in enumerate(fins):
    x1,y1 = np.loadtxt(f, delimiter=',',unpack=True)
    data.append( [f'poly-{n+1}',x1,y1] )

cmap = plt.get_cmap('jet')
colormap = cmap(np.linspace(0, 1, len(fins)))
color = [ colormap[c] for c in range(len(fins)) ] #'skyblue'
color = ['tab:orange','tab:blue']

fig, axs = plt.subplots(1,1,figsize=(3,1.5),tight_layout=True,)#dpi=100)

labels = ['Randomly placed GO','Vertically placed GO']
for n,d in enumerate(data):
    d[2] = d[2] / np.trapz( d[2], x=d[1] ) ## normalize by density

    #axs.plot(d[1],d[2], color=color[n], lw=1, ls='', marker='o', ms=2, label=d[0])
    y1 = denoise2(d[2])
    axs.plot(d[1],y1, color=color[n], lw=2, ls='-', marker='', ms=3, label=labels[n])

    axs.set_xlabel('GO orientation $(\degree)$',fontsize=10) ## input X name
    axs.set_ylabel('Distribution \ndensity',fontsize=10) ## input Y name
    axs.set_yticklabels([])
    axs.tick_params(direction='in',labelsize='8' )
    axs.set_ylim([0,0.04])
    axs.set_xlim([0,90])
    axs.set_xticks( range(0,91,15) )
    #axs.set_yticks( range(0,0.03,0.01) )


#plt.subplots_adjust(wspace=0,bottom=0.15)#, hspace=0.0)
axs.legend( loc='upper right', fontsize=9, frameon=False,)#  bbox_to_anchor=(0.99, 0.5))
plt.savefig('angle_go_1d.png', dpi=800)


exit()

