import mdtraj as md
import os
import matplotlib.pyplot as plt
from matplotlib import cm, colors
import sys
import numpy as np


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

    # BDO and water, and PO
    #bdo_atom = [at.index for at in topo.atoms if at.residue.name=='BDO' and at.element.symbol!='H' ]
    #wat_atom = [at.index for at in topo.atoms if at.residue.name=='HOH' and at.element.symbol!='H' ]
    bdo_atom = [at.index for at in topo.atoms if at.residue.name=='BDO' and at.name in ['O02','O05','H0F','H0A'] ]
    wat_atom = [at.index for at in topo.atoms if at.residue.name=='HOH' ]#and at.element.symbol=='O']

    po_N = [at.index for at in topo.atoms if at.residue.name==PO_name and at.element.symbol not in ['H','C']]
    go_C = [at.index for at in topo.atoms if at.residue.name=='GOF' and at.element.symbol not in ['H','C'] ]

    index_atom = np.array( po_N )
    index_odd = index_atom[::2] # first N
    index_even = index_atom[1::2] # second N

    """
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
    po_atom = np.transpose( [index_odd, index_even, c_ngb_atom] )  ## All pairs: N,N,C, 
    """

    # pairs for gr
    #pairs = [ [bdo_atom,c_ngb_atom], [wat_atom,c_ngb_atom] ]
    #pairs = [ [bdo_atom,index_odd], [wat_atom,index_odd] ]
    #pairs = [ [bdo_atom,index_even], [wat_atom,index_even] ]
    pairs = [ [wat_atom,po_N], [bdo_atom,po_N] ]
    all_gr = []
    for pair in pairs:
        pair = topo.select_pairs(selection1=pair[0], selection2=pair[1])
        print('Start gr calculation:')
        dr, gr= md.compute_rdf(traj, pair, r_range=(0.0, 1.0), bin_width=0.01)
        all_gr.append( [dr, gr] )

    if len(go_C)>0:
        pairs = [ [wat_atom,go_C], [bdo_atom,go_C] ]
        for pair in pairs:
            print('Start gr calculation:')
            pair = topo.select_pairs(selection1=pair[0], selection2=pair[1])
            dr, gr= md.compute_rdf(traj, pair, r_range=(0.0, 1.0), bin_width=0.01)
            all_gr.append( [dr, gr] )

    return all_gr


path_pbi = '/qfs/projects/sepcon/difan/sepcon-pbi/work-gmx/jobs-poly/'

gros = [
path_pbi+ 'job-msd01-poly1/nvt5.gro',
'job-box-poly1-5layeredGO/md-npt-03.gro',
'job-box-poly1-5randomGO/md-npt-03.gro',
#path_pbi+ 'job-msd01-poly3/nvt5.gro',
#'job-box-poly3-5layeredGO/md-npt-03.xtc',
#'job-box-poly3-5randomGO/md-npt-03.xtc',
path_pbi+ 'job-msd01-poly4/nvt5.gro',
'job-box-poly4-5layeredGO/md-npt-03.gro',
'job-box-poly4-5randomGO/md-npt-03.gro',
]
xtcs = [
path_pbi+ 'job-msd01-poly1/nvt8.xtc',
'job-box-poly1-5layeredGO/md-npt-03.xtc',
'job-box-poly1-5randomGO/md-npt-03.xtc',
#path_pbi+ 'job-msd01-poly3/nvt8.xtc',
#'job-box-poly3-5layeredGO/md-npt-03.xtc',
#'job-box-poly3-5randomGO/md-npt-03.xtc',
path_pbi+ 'job-msd01-poly4/nvt8.xtc',
'job-box-poly4-5layeredGO/md-npt-03.xtc',
'job-box-poly4-5randomGO/md-npt-03.xtc',
]

gros=gros[:3]
xtcs=xtcs[:3]   
#gros=gros[3:6]
#xtcs=xtcs[3:6]   
print( gros, xtcs )

# shape = gros, 4, 2
data = [ mdtj_gr(fxtc, fgro) for fgro,fxtc in zip(gros,xtcs) ] 

#label = [ i.split('/')[0][-5:] for i in gros ]
label = ['WAT-Poly','BDO-Poly','WAT-GO','BDO-GO' ]

cmap = plt.get_cmap('jet')
colormap = cmap(np.linspace(0, 1, len(gros)))
color = [ colormap[c] for c in range(len(gros)) ] #'skyblue'
color = ['r','g','b','m',]#'c','y','gray']
#color = ['r','seagreen','skyblue','m'] #+['b','gold','skyblue']+['m','gold','skyblue']

#n_plot = np.amax([len(d) for d in data ])
n_plot = len(data)
fig, axs = plt.subplots(1,1,figsize=(2,3.5), tight_layout=True)
axs.tick_params(direction='in',labelsize='8')
axs.set_xlim(left=0.1,right=0.7)
axs.set_ylim([-0.1, 6])
axs.set_yticklabels([])
#axs.set_xticklabels([])
axs.set_ylabel('g(r)',fontsize=10) ## input X name
axs.set_xlabel('r (nm)',fontsize=10) ## input Y name

y_shift = np.cumsum([0, 2, 2])
for i in range( len(data) ):  # for each poly system
    for j in range( len(data[i]) ): # for each pair
        xy = data[i][j]
        #axss[j].plot( xy[0],xy[1],linewidth=0.5,linestyle='-',marker='',markersize=1,color=color[i],label=label[i])
        if i==len(data)-1:
            axs.plot( xy[0],xy[1]+y_shift[i],linewidth=1,linestyle='-',marker='',markersize=1,color=color[j],label=label[j])
        else:
            axs.plot( xy[0],xy[1]+y_shift[i],linewidth=1,linestyle='-',marker='',markersize=1,color=color[j])
#axs.legend( loc='upper left', fontsize=10, frameon=False, ncol=2, columnspacing=2, labelspacing=0.25,)#  bbox_to_anchor=(0.01, 0.15))

#plt.subplots_adjust(wspace=0,hspace=0.0, bottom=0.15)

figout = 'fig-gr'+'.png'
plt.savefig(figout, format='png', dpi=800)
