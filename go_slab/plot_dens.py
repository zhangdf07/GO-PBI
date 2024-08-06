import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def read_xvg(fname):
    with open(fname, 'r') as ff1:
        lines = ff1.readlines()
        head = [ l for l in lines if l[0] in ["#","@"] ]
        body = [ l.split() for l in lines if l not in head ]
        data = [ [] for n in range(len(body[0])) ]
        for b in body:
            for n in range(len(data)):
                data[n].append( float(b[n]) )
    return head,data



fins = [
#'job-010-solvent-large',
'job-012-at-gr',
'job-013-groh-01',
'job-021-groh-02',
'job-014-groh-03',
'job-022-groh-04',
'job-015-groh-05',
'job-023-groh-06',
'job-024-groh-08',
'job-016-gro-01',
'job-025-gro-02',
'job-017-gro-03',
'job-026-gro-04',
'job-018-gro-05',
'job-027-gro-06',
'job-028-gro-08',
]


data = [
'dens-bdo.xvg',
'dens-wat.xvg',
]

allread = [[ read_xvg(f+"/"+d) for d in data ] for f in fins  ]

if False:
  colr = ['r','b']
  mark = ['v','o']
  fg = [0 for c in colr]
  lg = ['BDO','Water']
  nfin = len(fins)
  fig, axs = plt.subplots(1,nfin,figsize=(3*nfin,3),tight_layout=True)
  for i in range(nfin):
    axs[i].set_ylabel('Number density (1/nm$^3$)',fontsize=10) ## input Y name
    axs[i].set_xlabel('Z (nm)',fontsize=10) ## input X name
    axs[i].set_xlim(left=0,right=10)
    axs[i].tick_params(direction='in',labelsize='10')
    axs[i].set_ylim([0, 1200])
    for j in range(2):
      da = allread[i][j][1]
      if i in [5,6,7,8]:
        da[0] = np.array(da[0]) + 2
      fg[j],=axs[i].plot(da[0],da[1],linewidth=0.5,linestyle='-',marker=mark[j],markersize=1,color=colr[j])
    axs[i].legend(fg,lg,edgecolor='inherit',fontsize='10',loc=0,ncol=1,labelspacing=0.2,columnspacing=0.3,borderpad=0.2,handletextpad=0.3)

  figout = 'fig_density'+'.png'
  plt.savefig(figout, format='png', dpi=200)

if True:
  nfin = len(fins)
#  colr = ['r','b','g','m']
  colr = [ cm.tab20(n) for n in np.linspace(0,1,nfin) ]
  mk = ['o']+['v']*7+['s']*7
  fg = [0 for c in colr]
  lg = ['Gra-0%']
  lg += ['OH-10%','OH-20%','OH-30%','OH-40%','OH-50%','OH-60%','OH-80%']
  lg += ['O-10%','O-20%','O-30%','O-40%','O-50%','O-60%','O-80%']

  fig, axs = plt.subplots(1,2,figsize=(8,4),tight_layout=True)
  for n in range(2):
    axs[n].set_ylabel('Number density (1/nm$^3$)',fontsize=10) ## input Y name
    axs[n].set_xlabel('Z (nm)',fontsize=10) ## input X name
    axs[n].set_xlim(left=0,right=10)
    axs[n].set_ylim([0, 1200])
    axs[n].tick_params(direction='in',labelsize='10')

  d1 = [0] + list(range(1,8))
  d2 = [0] + list(range(8,nfin))

  #for n in range(nfin):
  for n in d1:
    for j in range(2):
      da = allread[n][j][1]
      fg[n],=axs[j].plot(da[0],da[1],linewidth=0.5,linestyle='-',marker=mk[n],markersize=1,color=colr[n])
  for n in range(2):
    axs[n].legend(fg,lg,edgecolor='inherit',fontsize='8',loc=1,ncol=1,labelspacing=0.2,columnspacing=0.3,borderpad=0.2,handletextpad=0.3)

  figout = 'fig_density'+'.png'
  plt.savefig(figout, format='png', dpi=200)

exit()
