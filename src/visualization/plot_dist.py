import os,sys, math, numpy as np, itertools
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from pylab import *
import src.utilities as utils

config = utils.read_config()
mpl.rcParams.update(mpl.rcParamsDefault)  # VS Code plots not black
plt.style.use(config['viz'])

infile='dist.dat'			  #First input file
outname='dist'	  		          #Name output files will take
xlbl='Amino Acid Number'
ylbl='Amino Acid Number'
ttl=''
maxc=16

mi=[]
mj=[]
ol=[]
i=-1

#############################################################################
# Read arguments from terminal, and assign input files and a name that all output files will contain. 
#############################################################################
for x in range(1,len(sys.argv)):
	if sys.argv[x] == '-i':
		infile = sys.argv[x+1]

	if sys.argv[x] == '-out':
		outname = sys.argv[x+1]
		
	if sys.argv[x]=='-xlabel':
		xlbl = sys.argv[x+1]
		
	if sys.argv[x]=='-ylabel':
		ylbl = sys.argv[x+1]
	
	if sys.argv[x]=='-title':
		ttl = sys.argv[x+1]
	
	if sys.argv[x]=='-val':
		maxc = sys.argv[x+1]
	
	if sys.argv[x]=='-help':
		print('\n\nProgram to plot overlap data...\n\nOPTIONS:\n'\
		'-i = Name of input file (Default=overlap.dat)\n'\
		'-xlabel = Label for x axis (Default=mode i)\n'\
		'-ylabel = Label for y axis (Default=mode j)\n'\
		'-title = Title for plot\n')
		exit()
		
inlines=open(infile,'r').readlines()

if inlines[-1]=='\n':
	inlines[-1:]=[]

i=i+1
mi.append([])
mj.append([])
ol.append([])

for line in inlines:
	if line=='\n':
		i=i+1
		mi.append([])
		mj.append([])
		ol.append([])
		
	else:
		mi[i].append(int(line.split()[0]))
		mj[i].append(int(line.split()[1]))
		ol[i].append(float(line.split()[2]))
		
mi=np.array(mi)
mj=np.array(mj)
ol=np.array(ol)

maxv = mi.max()
for x in range(1,len(sys.argv)):
	if sys.argv[x] == '-max':
		maxv = float(sys.argv[x+1])

fig=plt.figure(1, figsize=(11,8))
ax=fig.add_subplot(111)
cmain=ax.pcolor(mi,mj,ol,vmin=0, vmax=maxc,cmap=plt.cm.gist_yarg_r)
ax.set_title(ttl)

ax.set_xlabel(xlbl)
ax.set_xlim(mi.min(), maxv)

ax.set_ylabel(ylbl)
ax.set_ylim(maxv, mj.min())

cbar=fig.colorbar(cmain,aspect=10,ticks=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30])

fig.text(.85, .95, 'Distance / $\AA{}$', horizontalalignment='center')


# plt.rcParams.update({'font.size': 22})


plt.savefig(outname+'.png',format='png')
plt.show()

print('DONE')
