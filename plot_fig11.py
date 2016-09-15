import os.path

ldat = []
lsim = ['2014nov23_AMPAStimBatch_0_loadedrxd', '2014nov23_AMPAStimBatch_150_loadedrxd']


for sim in lsim:
		if os.path.isfile('./data/'+sim+'_.npz'):
				print 'found file:', sim
				ldat.append(myloaddata(sim))
if len(ldat)>0:
	lsec,lopt,lval = AMPAStimBatchParams()
	AMPAStimBatchDraw_fig11(dat=ldat)#,lidx=[0,1])
else: print "I couldn't open simulation data file(s). Make sure that you ran at least one simulation before trying to plot."





