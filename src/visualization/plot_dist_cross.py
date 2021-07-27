import os,sys, math, numpy as np, itertools
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from pylab import *
from numpy import ma
import src.utilities as utils

config = utils.read_config()
mpl.rcParams.update(mpl.rcParamsDefault)  # VS Code plots not black
plt.style.use(config['viz'])

infile='dist.dat'			  #First input file
infile2='crosscor.dat'                    #Second input file
outname='graph'	  		          #Name output files will take
xlbl='Amino Acid Number'
ylbl='Amino Acid Number'
ttl=''

mi=[]
mi2=[]
mj=[]
mj2=[]
ol=[]
ol2=[]
i=-1
j=-1

#############################################################################
# Read arguments from terminal, and assign input files and a name that all output files will contain. 
#############################################################################
for x in range(1,len(sys.argv)):
	if sys.argv[x] == '-dist':
		infile = sys.argv[x+1]

	if sys.argv[x] == '-cross':
		infile2 = sys.argv[x+1]
		
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

#------------------------------------------------------------------

inlines=open(infile2,'r').readlines()

if inlines[-1]=='\n':
	inlines[-1:]=[]

j=j+1
mi2.append([])
mj2.append([])
ol2.append([])

for line in inlines:
	if line=='\n':
		j=j+1
		mi2.append([])
		mj2.append([])
		ol2.append([])
		
	else:
		mi2[j].append(int(line.split()[0]))
		mj2[j].append(int(line.split()[1]))
		ol2[j].append(float(line.split()[2]))
		
mi2=np.array(mi2)
mj2=np.array(mj2)
ol2=np.array(ol2)


#---------------------------------------------------------------------

Zm = ma.masked_where(mi2 < mj2, ol2)

fig=plt.figure(1, figsize=(13,9))
ax=fig.add_subplot(111,autoscale_on=False)
cmain=ax.pcolor(mi,mj,ol,vmin=0, vmax=16,cmap=plt.cm.gist_yarg_r)
csec=ax.pcolor(mi2,mj2,Zm,vmin=-1, vmax=1,cmap=plt.cm.RdBu_r)
ax.set_aspect(1)

ax.set_title(ttl)

ax.set_xlabel(xlbl)
ax.set_xlim(mi.min(), mi.max())

ax.set_ylabel(ylbl)
ax.set_ylim(mj.max(), mj.min())

cbar=fig.colorbar(cmain,aspect=10,shrink=0.9,ticks=[0,2,4,6,8,10,12,14,16])
cbar2=fig.colorbar(csec,aspect=10,shrink=0.9,ticks=[-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1])

gca()

for item in range(mi.min(), mi.max()):
	if item not in mi:
		gca().add_patch(Rectangle((item,mj.min()),1,mj.max()-mj.min(),color='darkseagreen'))
		gca().add_patch(Rectangle((mi.min(),item),mi.max()-mi.min(),1,color='darkseagreen'))

fig.text(.81, .9, 'Distance \n / $\AA{}$', horizontalalignment='center')
fig.text(.70, .9, 'Cross \n correlation', horizontalalignment='center')

# plt.rcParams.update({'font.size': 22})

plt.savefig(outname+'.png',format='png')
plt.show()

#-------------------------------------------------------------------

print('DONE')

