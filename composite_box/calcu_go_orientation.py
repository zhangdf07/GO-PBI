import mdtraj as md
import os
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import sys
import numpy as np


def unit_vec(vec):
    return vec/np.linalg.norm(vec)

def get_angle_of_vec(vec1,vec2):
    return np.arccos(  np.dot( unit_vec(vec1), unit_vec(vec2) )  )


def mdtj_gr(fxtc,fgro):
    # Read trajectory 
    print('Start reading:')
    traj = md.load(fxtc, top=fgro, stride=10)
    print('Old',fxtc, traj)
    topo = traj.topology
    #traj = traj[-100:]
    print('New',fxtc, traj)

    natom = topo.n_atoms
    nresd = topo.n_residues

    # Build specific atoms 
    atom = [at for at in topo.atoms if at.residue.name=='GOF' and at.element.symbol=='C'] 
    go_flakes = {}

    for at in atom:
        k = at.residue.index
        if k not in list( go_flakes.keys() ):
            go_flakes[k] = [at.index]
        else:
            go_flakes[k].append( at.index )
    
    go_flakes = [ v[:229 ] for v in go_flakes.values() ]  ## All the C in the GO orginal

    # Get neighbors
    cutoff = 0.2 
    bins1 = np.linspace(0,180,num=181)
    angle_distribution1 = []
    for go in go_flakes: # each GO
        print( 'Start GO' )
        # ngb of this go's atoms in all frames
        ngb = []
        for at in go:
            ngb += [ b for b in  md.compute_neighbors(traj[0], cutoff, [at], haystack_indices=go) if len(b)==3 ]
        ngb = np.array(ngb)
        # Get vectors: 
        dist1 = md.compute_displacements( traj, ngb[:,[0,1]] ) # periodic=True, i frame, j pair, 3
        dist2 = md.compute_displacements( traj, ngb[:,[0,2]] ) # periodic=True, i frame, j pair, 3
        for frame_indx in range( len(traj) ):  # loop all frames
            vec2 = np.cross( dist1[frame_indx], dist2[frame_indx] )
            vec2 = vec2/np.reshape( np.linalg.norm(vec2, axis=1), (len(vec2),-1) ) # The vectors in each frame
            vec3 = np.copy(vec2)
            vec3[:,0]=0

            angle = []
            for v2,v3 in zip( vec2,vec3 ):
                prod1= np.dot(v2,v3)
                if prod1>1:
                    prod1 = 1
                elif prod1<-1:
                    prod1 = -1
                angle.append( np.arccos( prod1 ) )
            angle = np.degrees( angle )

            mask = angle>90
            angle = np.abs( angle - mask*180 )

            hist1, edges1 = np.histogram( angle, bins=bins1 )
            angle_distribution1.append( hist1 )

    #print( angle_distribution1 )
    angle_distribution1 = np.mean( angle_distribution1, axis=0)
    edges1 = (edges1[:-1]+edges1[1:])/2

    return edges1, angle_distribution1


gros = [
'job-box-poly1-5layeredGO/md-npt-03.gro',
'job-box-poly1-5randomGO/md-npt-03.gro',
#'job-box-poly3-5layeredGO/md-npt-03.xtc',
#'job-box-poly3-5randomGO/md-npt-03.xtc',
'job-box-poly4-5layeredGO/md-npt-03.gro',
'job-box-poly4-5randomGO/md-npt-03.gro',
]
xtcs = [
'job-box-poly1-5layeredGO/md-npt-03.xtc',
'job-box-poly1-5randomGO/md-npt-03.xtc',
#'job-box-poly3-5layeredGO/md-npt-03.xtc',
#'job-box-poly3-5randomGO/md-npt-03.xtc',
'job-box-poly4-5layeredGO/md-npt-03.xtc',
'job-box-poly4-5randomGO/md-npt-03.xtc',
]

gros=gros[:2]
xtcs=xtcs[:2]

label = [ i.split('/')[0][8:] for i in gros ]

cmap = plt.get_cmap('jet')
colormap = cmap(np.linspace(0, 1, len(gros)))
color = [ colormap[c] for c in range(len(gros)) ] #'skyblue'

## The 1d distribution plot
fig, axs = plt.subplots(1,1,figsize=(3,3),tight_layout=True,dpi=100)

for n,(fgro,fxtc) in enumerate(zip(gros,xtcs)):
    #fxtc = fgro
    x1,y1 = mdtj_gr(fxtc, fgro)
    axs.plot(x1,y1, color=color[n], lw=1, ls='-', marker='', ms=3, label=label[n])

    output = np.transpose([x1,y1])
    np.savetxt(f"go_orientation-{label[n]}.csv", output, delimiter=",")

plt.legend( loc='upper right', fontsize=10, frameon=False,)#  bbox_to_anchor=(0.99, 0.5))
plt.savefig('angle_go_1d.png', dpi=200)


exit()


# The 2d distribution plot
fig, axs = plt.subplots(1,len(xy_2d),figsize=(3*len(xy_2d),3),tight_layout=True,dpi=100)
for n,xy in enumerate( xy_2d.values() ):
    #axs[n].plot(xy[0],xy[1], color=color[n], lw=1, ls='', marker='o', ms=3, label=label[n], alpha=0.6)
    hist = axs[n].hist2d(xy[0],xy[1], bins=180, )#norm=colors.LogNorm() )

plt.savefig('angle_distribution_2d.png', dpi=200)

#import pickle 
#with open('saved_xy_2d.pkl', 'wb') as f:
#    pickle.dump(xy_2d, f)

