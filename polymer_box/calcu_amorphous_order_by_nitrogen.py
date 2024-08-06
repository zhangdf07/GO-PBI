import mdtraj as md
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import numpy as np


def unit_vec(vec):
    return vec/np.linalg.norm(vec)

def get_angle_of_vec(vec1,vec2):
    return np.arccos(  np.dot( unit_vec(vec1), unit_vec(vec2) )  )


def mdtj_gr(fxtc,fgro):
    print('Start reading:')
    traj = md.load(fxtc, top=fgro, stride=2)
    print('Old',fxtc, traj)
    topo = traj.topology
    traj = traj[-300:]
    print('New',fxtc, traj)

    natom = topo.n_atoms
    nresd = topo.n_residues

    print('Building pairs')

    PO_name = set( [at.residue.name for at in topo.atoms if 'PO' in at.residue.name] )
    PO_name = list(PO_name)[0]

    atom = [at for at in topo.atoms if at.residue.name==PO_name and at.element.name=='nitrogen']
    #print( atom )

    """ 
    # Only when we consider residuces
    residue =  np.sort( np.array( list(set([at.residue for at in atom ])), dtype=str ) )

    dict_residue = { a:[] for a in residue }
    for at in atom:
        k = np.array([at.residue], dtype=str)[0]
        dict_residue[ k ].append(at.index)

    #po_atom = np.transpose([ v for v in dict_residue.values() ])  ## loop all atoms in eempa, then loop all eempa, each value is atom index
    po_atom = np.array([ v for v in dict_residue.values() ])  ## loop all molecules, then loop all atoms, each value is atom index
    """
    index_atom = np.array( [ at.index for at in atom ] )
    index_odd = index_atom[::2]
    xyz_odd = traj.xyz[:,index_odd,:] # i frame, j atom, 3 for X Y Z

    # calculate local neighbor
    # index_odd is the starting point of all vectors
    cutoff = 1 # nm
    ngb_list = {}
    for i in index_odd:
        others = index_odd[ np.where(index_odd!=i) ]
        others = np.transpose( [np.full(len(others),i), others] )  # [ [i,j1],[i,j2]... ]
        dist = md.compute_distances(traj,others)  ## i frames, k pairs of i-j

        ngb_i = [] ## the ngb of i in each frame
        for d in dist: # every frame:
            mask = d<cutoff
            ngb =  others[mask][:,1]
            ngb_i.append(ngb)

        ngb_list[i] = ngb_i

    ## The we go to vector
    po_atom = np.reshape( index_atom, (-1,2) )  ## The pairs of N: j pairs, 2
    dist = md.compute_displacements(traj, po_atom) # periodic=True, i frame, j pair, 3

    bins = np.linspace(0,90,num=91) ## 0,1,2...90 if num=91
    angle_distribution = []
    for frame_indx, displacement in enumerate(dist):  # loop frames
        #xyz = xyz + dsp/2  # get the center of vec, i.e. mid point of each N in the pair. This is the starting point XYZ + 1/2 displacement

        # change these vectors to unit vectors (dsp)
        vec_length = np.reshape( np.linalg.norm(displacement, axis=1), (len(displacement),-1) )
        dsp = displacement / vec_length
        """
        # Calculate angles for all vector
        angles_in_frame = []
        for i in range( 0, len(dsp)-1 ):
            for j in range( i+1, len(dsp) ):
                prod = np.dot(dsp[i],dsp[j])
                if prod>1:
                    prod = 1
                elif prod<-1:
                    prod = -1
                angles_in_frame.append( np.arccos( prod ) )
        """
        # Calculate angles for only ngbs
        angles_in_frame = []
        for i in range(len(dsp)):
            ngb = ngb_list[ index_odd[i] ][frame_indx] # The ngb of vec starting point (index_odd[i]) in this frame index
            for ngb_j in ngb:
                j = np.argwhere(index_odd==ngb_j).flatten()[0]
                prod = np.dot(dsp[i],dsp[j])
                if prod>1:
                    prod = 1
                elif prod<-1:
                    prod = -1
                angles_in_frame.append( np.arccos( prod ) )
         
        # Convert to degrees       
        angles_in_frame = np.degrees(angles_in_frame)

        mask = angles_in_frame>90
        angles_in_frame = np.abs( angles_in_frame - mask*180 ) ## angle>90, then use 180-angle so that all angles are 0<a<90
        #print( angles_in_frame, angles_in_frame.shape  )

        # Bin the angle distribution
        hist, edges = np.histogram( angles_in_frame, bins=bins )
        angle_distribution.append( hist )

    angle_distribution = np.mean( angle_distribution, axis=0)
    angle_distribution = angle_distribution/np.linalg.norm(angle_distribution)
    edges = (edges[:-1]+edges[1:])/2

    return edges, angle_distribution

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

label = [ i.split('/')[0][-5:] for i in gros ]

cmap = plt.get_cmap('jet')
colormap = cmap(np.linspace(0, 1, len(gros)))
color = [ colormap[c] for c in range(len(gros)) ] #'skyblue'

fig, axs = plt.subplots(1,1,figsize=(3,3),tight_layout=True,dpi=100)

base_line = None
for n,(fgro,fxtc) in enumerate(zip(gros,xtcs)):
    #fxtc = fgro
    x,y = mdtj_gr(fxtc, fgro)
    if base_line is None:
        base_line = np.copy(y)
    #axs.bar(x,y, color=color[n], width=0.5, label=label[n])
    #axs.plot(x,y-base_line, color=color[n], lw=1, ls='-', marker='', ms=3, label=label[n])
    axs.plot(x,y, color=color[n], lw=1, ls='-', marker='', ms=3, label=label[n])

    output = np.transpose([x,y])
    np.savetxt(f"angle-1nm-{n+1}.csv", output, delimiter=",")

plt.legend( loc='upper left', fontsize=10, frameon=False,)#  bbox_to_anchor=(0.99, 0.5))

plt.savefig('angle_distribution.png', dpi=200)

