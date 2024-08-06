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
    traj = md.load(fxtc, top=fgro, stride=2)
    print('Old',fxtc, traj)
    topo = traj.topology
    traj = traj[-300:]
    print('New',fxtc, traj)

    natom = topo.n_atoms
    nresd = topo.n_residues

    # Build specific atoms 
    print('Building pairs')
    PO_name = set( [at.residue.name for at in topo.atoms if 'PO' in at.residue.name] )
    PO_name = list(PO_name)[0]

    atom = [at for at in topo.atoms if at.residue.name==PO_name and at.element.name=='nitrogen'] # N in polymer
    index_atom = np.array( [ at.index for at in atom ] )
    index_odd = index_atom[::2] # first N
    index_even = index_atom[1::2] # second N
    resid = [ topo.atom(i).residue.index for i in index_odd ]

    # Find the C ngb of N using the first frame of traj
    cutoff = 0.2 
    c_ngb_atom = [] 
    for p1,p2 in zip(index_odd,index_even):
        ngb0 = md.compute_neighbors(traj, cutoff, [ p1 ]) # N1's ngb
        ngb1 = md.compute_neighbors(traj, cutoff, [ p2 ]) # N2's ngb
        # only cosider the first frame
        c_ngb = np.intersect1d( ngb0[0], ngb1[0] ) # Commone ngb of N1 and N2 is the C we need
        c_ngb = np.array([ c for c in c_ngb if topo.atom(c).element.name=='carbon' ]) # Sometimes the terminal H is also counted
        if len(c_ngb)!=1: # There should be only one C
            print( c_ngb , 'is wrong in', fxtc)
            exit()
        else:
            c_ngb_atom.append( c_ngb.item() )
    c_ngb_atom = np.array( c_ngb_atom ) # third C
    po_atom = np.transpose( [index_odd, index_even, c_ngb_atom, resid] )  ## All pairs: N,N,C, resid

    # Get vectors: N->C
    dist1 = md.compute_displacements( traj, po_atom[:,[2,0]] ) # periodic=True, i frame, j pair, 3
    dist2 = md.compute_displacements( traj, po_atom[:,[2,1]] ) # periodic=True, i frame, j pair, 3

    # Get the local neighbor list using C atom
    cutoff = 0.5 # 1 # nm
    ngb_list = {}  # key is index of C. value is array: i frame, x ngbs
    for i in c_ngb_atom:
        """
        others = c_ngb_atom[ np.where(c_ngb_atom!=i) ]  # The ngb of i
        others = np.transpose( [np.full(len(others),i), others] )  # [ [i,j1],[i,j2]... ]
        dist = md.compute_distances(traj,others)  ## i frames, j pairs of i-j
        ngb_i = [ others[d<cutoff][:,1] for d in dist ]  ## The ngb of i in all frames
        """
        ngb_i = md.compute_neighbors(traj, cutoff, [i], haystack_indices=c_ngb_atom )
        ngb_list[i] = ngb_i  

    print('Computing angles')
    # Calculate angles between two vectors: angle 1 is vec between dist1+dist2. angle 2 is vec between normal vec (plane orientation)
    bins1 = np.linspace(0,180,num=181) ## 0,1,2...180 if num=181
    bins2 = np.linspace(0,180,num=181) ## 0,1,2...180 if num=181
    #bins2 = np.linspace(0,90,num=91) ## 0,1,2...90 if num=91
    angle_distribution1, angle_distribution2 = [],[]
    xy_distribution = [ [], [] ] 
    for frame_indx in range( len(traj) ):  # loop all frames
        # For vectors we want: 1=planar and 2=vertical, into unit vector. Size is i frame* j pair
        vec1 = dist1[frame_indx] + dist2[frame_indx]  
        vec2 = np.cross( dist1[frame_indx], dist2[frame_indx] ) 
        vec1 = vec1/np.reshape( np.linalg.norm(vec1, axis=1), (len(vec1),-1) )
        vec2 = vec2/np.reshape( np.linalg.norm(vec2, axis=1), (len(vec2),-1) )

        # Only calculate ngb's angles
        angles_in_frame_1,angles_in_frame_2 = [],[]  ## in this frame
        for i in range(len(po_atom)): ## loop all the NNC pairs-
            ngb = ngb_list[ c_ngb_atom[i] ][frame_indx] # the ngb of c_ngb_atom[i] in this frame
            #ngb = [ c for c in c_ngb_atom if c!=c_ngb_atom[i] ]  # all c
            for ngb_j in ngb:
                j = np.argwhere(c_ngb_atom==ngb_j).item() # Get the index of this ngb so that we know where to find its vec12
                if  resid[i] != resid[j]:
                    prod1 = np.dot(vec1[i],vec1[j])
                    prod2 = np.dot(vec2[i],vec2[j])
                    if prod1>1:
                        prod1 = 1
                    elif prod1<-1:
                        prod1 = -1
                    if prod2>1:
                        prod2 = 1
                    elif prod2<-1:
                        prod2 = -1
                    angles_in_frame_1.append( np.arccos( prod1 ) )
                    angles_in_frame_2.append( np.arccos( prod2 ) )

        # Convert to degrees       
        angles_in_frame_1 = np.degrees(angles_in_frame_1)
        angles_in_frame_2 = np.degrees(angles_in_frame_2)
        #mask = angles_in_frame_2>90  # for planar orientation, angles are no larger than 90
        #angles_in_frame_2 = np.abs( angles_in_frame_2 - mask*180 ) ## if angle>90, then use 180-angle so that all angles are 0<a<90
        xy_distribution[0] += list(angles_in_frame_1)
        xy_distribution[1] += list(angles_in_frame_2)

        # Bin the angle distribution
        hist1, edges1 = np.histogram( angles_in_frame_1, bins=bins1 )
        hist2, edges2 = np.histogram( angles_in_frame_2, bins=bins2 )
        angle_distribution1.append( hist1 )
        angle_distribution2.append( hist2 )

    angle_distribution1 = np.mean( angle_distribution1, axis=0)
    angle_distribution2 = np.mean( angle_distribution2, axis=0)
    #angle_distribution1 = angle_distribution1/np.linalg.norm(angle_distribution1)
    #angle_distribution2 = angle_distribution2/np.linalg.norm(angle_distribution2)
    edges1 = (edges1[:-1]+edges1[1:])/2
    edges2 = (edges2[:-1]+edges2[1:])/2

    return edges1, angle_distribution1, edges2, angle_distribution2, xy_distribution

gros = [
 'job-msd01-poly1/nvt5.gro',
 'job-msd01-poly2/nvt5.gro',
 'job-msd01-poly3/nvt5.gro',
 'job-msd01-poly4/nvt5.gro',
 'job-msd01-poly5/nvt5.gro',
 'job-msd01-poly6/nvt5.gro',
 'job-msd01-poly7/nvt5.gro',
]
xtcs = [
 'job-msd01-poly1/nvt8.xtc',
 'job-msd01-poly2/nvt8.xtc',
 'job-msd01-poly3/nvt8.xtc',
 'job-msd01-poly4/nvt8.xtc',
 'job-msd01-poly5/nvt6.xtc',
 'job-msd01-poly6/nvt6.xtc',
 'job-msd01-poly7/nvt6.xtc',
]

#gros=gros[:4]
#xtcs=xtcs[:4]

label = [ i.split('/')[0][-5:] for i in gros ]

cmap = plt.get_cmap('jet')
colormap = cmap(np.linspace(0, 1, len(gros)))
color = [ colormap[c] for c in range(len(gros)) ] #'skyblue'

## The 1d distribution plot
fig, axs = plt.subplots(1,2,figsize=(7,3),tight_layout=True,dpi=100)

xy_2d = {}
base_line = None
for n,(fgro,fxtc) in enumerate(zip(gros,xtcs)):
    #fxtc = fgro
    x1,y1,x2,y2, xy = mdtj_gr(fxtc, fgro)
    xy_2d[n] = xy 
    axs[0].plot(x1,y1, color=color[n], lw=1, ls='-', marker='', ms=3, label=label[n])
    axs[1].plot(x2,y2, color=color[n], lw=1, ls='-', marker='', ms=3, label=label[n])

    output = np.transpose([x1,y1,x2,y2])
    np.savetxt(f"angle_1d-05nm-{n+1}.csv", output, delimiter=",")

    output = np.transpose(xy)
    np.savetxt(f"angle_2d-05nm-{n+1}.csv", output, delimiter=",")

plt.legend( loc='upper left', fontsize=10, frameon=False,)#  bbox_to_anchor=(0.99, 0.5))
plt.savefig('angle_distribution_1d.png', dpi=200)

# The 2d distribution plot
fig, axs = plt.subplots(1,len(xy_2d),figsize=(3*len(xy_2d),3),tight_layout=True,dpi=100)
for n,xy in enumerate( xy_2d.values() ):
    #axs[n].plot(xy[0],xy[1], color=color[n], lw=1, ls='', marker='o', ms=3, label=label[n], alpha=0.6)
    hist = axs[n].hist2d(xy[0],xy[1], bins=180, )#norm=colors.LogNorm() )

plt.savefig('angle_distribution_2d.png', dpi=200)

#import pickle 
#with open('saved_xy_2d.pkl', 'wb') as f:
#    pickle.dump(xy_2d, f)

