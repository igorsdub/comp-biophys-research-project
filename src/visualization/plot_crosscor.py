import os,sys, math, numpy as np, itertools
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from pylab import *
import src.utilities as utils

config = utils.read_config()
mpl.rcParams.update(mpl.rcParamsDefault)  # VS Code plots not black
plt.style.use(config['viz'])

infile='crosscor.dat'			  #First input file
outname='crosscor'	  		          #Name output files will take
xlbl='Amino Acid Number'
ylbl='Amino Acid Number'
ttl=''

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
		
	if sys.argv[x]=='-xlabel':
		xlbl = sys.argv[x+1]
		
	if sys.argv[x]=='-ylabel':
		ylbl = sys.argv[x+1]
	
	if sys.argv[x]=='-title':
		ttl = sys.argv[x+1]
		
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

fig=plt.figure(1, figsize=(11,8))
ax=fig.add_subplot(111)
cmain=ax.pcolor(mi,mj,ol,vmin=-1, vmax=1,cmap=plt.cm.RdBu_r)
ax.set_title(ttl)
ax.set_xlabel(xlbl)
#ax.set_xticks([100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360])
#ax.set_xticklabels(['\n $100$', '\n $120$', '\n $140$', '\n $160$', '\n $180$', '\n $200$', '\n $220$', '\n $240$', '\n $260$', '\n $280$', '\n $300$', '\n $320$', '\n $340$', '\n $360$'])
ax.set_xlim(mi.min(), mi.max())

ax.set_ylabel(ylbl)
#ax.set_ylim(log(0.25), log(4))
#ax.set_yticks([log(0.25), log(0.5), log(0.75), log(1), log(2), log(3), log(4)])
#ax.set_yticklabels(['$0.25$', '$0.5$', '$0.75$', '$1$', '$2$', '$3$', '$4$'])
ax.set_ylim(mj.max(), mj.min())

cbar=fig.colorbar(cmain,aspect=10,ticks=[-1,-0.75,-0.5,-0.25,0.0,0.25,0.5,0.75,1.0])
#cbar.ax.set_yticklabels(['$0.965$','$0.975$','$0.985$','$0.995$','$1.005$','$1.015$','$1.025$'])

gca()

for item in range(mi.min(), mi.max()):
	if item not in mi:
		gca().add_patch(Rectangle((item,mj.min()),1,mj.max()-mj.min(),color='blue'))
		gca().add_patch(Rectangle((mi.min(),item),mi.max()-mi.min(),1,color='blue'))


fig.text(.85, .95, 'Cross correlation', horizontalalignment='center')

# plt.rcParams.update({'font.size': 22})


plt.savefig(outname+'.png',format='png')
plt.show()

print('DONE')
