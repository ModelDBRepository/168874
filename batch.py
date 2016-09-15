from pylab import *
import sys,os,numpy,subprocess
from neuron import h,rxd
from math import ceil
h.load_file("stdrun.hoc") # creates cvode object - not really needed here
from vector import *
from nqs import *
from conf import *
from cawave import myloaddata, caCYT_init, wavespeed, wavedur, waveonset, wavenq, wavetime, wavedist
import multiprocessing
from Queue import Queue
import matplotlib.gridspec as gridspec
#from IPython.core.debugger import Tracer; debug_here = Tracer() #debugging using ipdb

ion() # interactive mode

# append line s to filepath fn
def appline (s,fn):
  fp = open(fn,"a"); fp.write(s + "\n"); fp.close()

# check that the batch dir exists
def checkdir (d):
  try:
    if not os.path.exists(d): os.mkdir(d)
    return True
  except:
    print "could not create directory :" + d
    return False

# make a list of the sims that have already had their output saved, can then
# pass it into batchRun to skip those sims
def getSkipList (whichParams):
  lsec,lopt,lval = whichParams()
  sidx,lskip = -1,[]
  for i in xrange(len(lopt[0])):
    if lopt[0][i] == 'simstr':
      sidx = i
      break
  if sidx == -1:
    print "no simstr found!"
    return None
  for i in xrange(len(lval)):
    if os.path.exists("./data/" + lval[i][sidx] + "_.npz"):
      lskip.append(i)
  return lskip

# run a batch using multiprocessing 
#  based on http://www.bryceboe.com/2011/01/28/the-python-multiprocessing-queue-and-large-objects/
def batchRun (whichParams,blog,skip=[],qsz=10,bdir="./batch"):
  if not checkdir(bdir): return False
  jobs = multiprocessing.Queue()
  lsec,lopt,lval = whichParams()

  def myworker (jobs):
    while True:
      scomm = jobs.get()
      if scomm == None: break
      print "worker starting : " , scomm
      os.system(scomm) #worker function, invoked in a process.

  for i in xrange(len(lsec)):
    if i in skip: continue
    cfgname = os.path.join(bdir, str(i) + ".cfg")
    writeconf(cfgname,sec=lsec[i],opt=lopt[i],val=lval[i])
    cmd = "python cawave.py " + cfgname
    appline(cmd,blog)
    jobs.put(cmd)
  workers = []
  for i in range(qsz):
    jobs.put(None)
    tmp = multiprocessing.Process(target=myworker, args=(jobs,))
    tmp.start()
    workers.append(tmp)
  for worker in workers: worker.join()
  return jobs.empty()

# pass in a function pointer that specifies parameters (whichParams) and run the sims
def shortRun (whichParams,skip=[]):
  lsec,lopt,lval = whichParams()
  for i in xrange(len(lsec)):    
    if i in skip: continue
    writeconf("tmp.cfg",sec=lsec[i],opt=lopt[i],val=lval[i])
    os.system("python cawave.py tmp.cfg")

# return results from sims specified by function pointer whichParams 
def shortRead (whichParams):
  lsec,lopt,lval = whichParams()
  dat,sidx = [],0
  for i in xrange(len(lopt[0])):
    if lopt[0][i] == 'simstr':
      sidx=i
  for l in lval: dat.append(myloaddata(l[sidx])) # assumes that l[0] has simulation string
  return dat

# convert params to NQS
def params2nq (whichParams):
  lsec,lopts,lval = whichParams()
  ncol,nrow = len(lopts[0]),len(lopts)
  nq = NQS(ncol)
  for i in xrange(ncol): nq.s[i].s = lopts[0][i]
  if nq.fi("simstr") != -1: nq.strdec("simstr")
  for i in xrange(nrow):    
    for j in xrange(ncol):
      print i,j
      if nq.s[j].s == "simstr":
        nq.appi(j,lval[i][j])
      else:
        nq.appi(j,double(lval[i][j]))
  return nq

# load/return data from the simstr entry of nq at specified row
def loadfromnq (nq,row):
  nq.tog("DB")
  if row < 0 or row >= nq.v[0].size():
    print "row", row, " out of bounds: ", nq.v[0].size()
  simstr = nq.get("simstr",row).s
  print "loading " , simstr , " data "
  return myloaddata(simstr)

# add a column of wavenq objects
def addwavenqcol (nq,thresh):  
  nq.tog("DB"); 
  if nq.fi("nqw") == -1:
    nq.resize("nqw")
    nq.odec("nqw"); 
    nq.pad()
  sz = int(nq.v[0].size())
  for row in xrange(sz):
    print "up to row ", row, " of " , sz
    s = nq.get("simstr",row).s
    dat = loadfromnq(nq,row)
    nqw = wavenq(dat["cytca"],thresh)
    nq.set("nqw",row,nqw)
    del dat

# add the wave properties to the nqs
def addwavepropcols (nq):
  nq.tog("DB")
  if nq.fi("nqw") == -1:
    print "error: make sure this nqs has nqw (wavenq) column!"
    return False
  if nq.fi("speed") == -1: 
    nq.resize("speed","dur","time","amp","onset","dist"); nq.pad()
  sz = int(nq.v[0].size())
  for row in xrange(sz):    
    dat = loadfromnq(nq,row)
    nqw = nq.get("nqw",row).o[0]
    nqw.verbose=0
    nq.getcol("speed").x[row]=wavespeed(data=dat,nqw=nqw)
    nq.getcol("dur").x[row]=wavedur(data=dat,nqw=nqw); 
    nq.getcol("time").x[row]=wavetime(data=dat,nqw=nqw)
    nq.getcol("amp").x[row]=amax(dat["cytca"])
    nq.getcol("onset").x[row]=waveonset(data=dat,nqw=nqw)
    nq.getcol("dist").x[row]=wavedist(data=dat,nqw=nqw)
    nqw.tog("DB")
    nqw.verbose=1
    del dat
  return True
    

# append to the lists
def NewParam (lsec,lopt,lval,sec,opt,val):
  lsec.append(sec); lopt.append(opt); lval.append(val)

#####################################################################################
#                        baseline figure
# params for baseline figure
def baseParams ():
  lsec,lopt,lval = [],[],[]
  # baseline
  lsec.append(["run","run","run","run"])
  lopt.append(["simstr","runit","saveout","tstop"])
  lval.append(["simBase_","1","1","6000"])
  return lsec, lopt, lval

def baseRun ():
  shortRun(baseParams)

def baseRead ():
  return shortRead(baseParams)

# draw baseline simulation
def baseDraw (dat=None,miny=250,maxy=750):
  if dat is None: dat = baseRead()[0]
  from cawave import wavenq,mydraw,data,recdt
  mint,maxt=2,6
  tx, ty = -0.15, 1.06
  figure(); 
  ax=subplot(1,2,1); imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,6,0,1000),origin='lower'); colorbar();
  xlim((mint,maxt)); ylim((miny,maxy)); title(r'$Ca^{2+}_{cyt}$ (mM)'); xlabel('Time (s)'); ylabel(r'Position ($\mu$m)');
  text(tx,ty,"a",fontsize=14,fontweight='bold',transform=ax.transAxes)
  ax=subplot(1,2,2); imshow(dat["erca"],vmin=0,vmax=0.011,aspect='auto',extent=(0,6,0,1000),origin='lower'); colorbar();
  xlim((mint,maxt)); ylim((miny,maxy)); title(r'$Ca^{2+}_{er}$ (mM)'); xlabel('Time (s)'); 
  text(tx,ty,"b",fontsize=14,fontweight='bold',transform=ax.transAxes)


#####################################################################################

#####################################################################################
#                    slide 18 -- hotspots == extra IP3R or extra ER ?
# params for slide 18 ( baseline, double IP3R, quadruple ER )
def hsTestParams ():
  lsec,lopt,lval = [],[],[]
  # baseline
  lsec.append(["run","run","run"])
  lopt.append(["simstr","runit","saveout"])
  lval.append(["hsTestBase_","1","1"])
  # double IP3R
  lsec.append(["run","run","run","set","set"])
  lopt.append(["simstr","runit","saveout","ip3_origin","boost_every"])
  lval.append(["hsTest2IP3R_","1","1",str(2.0*12040.0),"50.0"])
  # quadruple ER
  lsec.append(["run","run","run","set","set"])
  lopt.append(["simstr","runit","saveout","boost_every","er_scale"])
  lval.append(["hsTest4ER_","1","1","50.0","4.0"])
  return lsec, lopt, lval

# run sims for slide 18 and save output ( baseline, double IP3R, quadruple ER )
def hsTestRun (skip=[]):
  shortRun(hsTestParams,skip)

# read results from slide 18 sims
def hsTestRead ():
  return shortRead(hsTestParams)

# draw results from slide 18
def hsTestDraw (dat=None):
  from cawave import wavenq,mydraw,data,recdt
  if dat is None: dat = hsTestRead()
  tx, ty = -0.15, 1.06
  sca = r'$Ca^{2+}_{cyt}$ (mM)'
  ltitles = ["Baseline: "+sca, r"$2X$ IP3R: "+sca, r"$4X$ ER: "+sca];  txl = ["a", "b", "c"];
  for i in xrange(len(ltitles)):
    nq = wavenq(dat[i]["cytca"]); nq.select("overlap",0,"up",1,"y",">=",500); ld=nqs2pyd(nq); nqsdel(nq);
    ax=subplot(2,len(ltitles),i+1); 
    imshow(dat[i]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower');
    xlim((0,12)); ylim((400,600));
    if i == 0: ylabel('Position (um)');
    # if i == len(ltitles)-1: colorbar()
    title(ltitles[i]); 
    text(tx,ty,txl[i],fontsize=14,fontweight='bold',transform=ax.transAxes)
    subplot(2,len(ltitles),i+4);
    plot(ld['startt']/1e3,ld['durt']/1e3,'k');  plot(ld['startt']/1e3,ld['durt']/1e3,'ko');
    # plot(ld['startt']/1e3,ld['speed']*500,'r');  plot(ld['startt']/1e3,ld['speed']*500,'ro');
    xlabel('Time(s)'); 
    if i==0: ylabel('Duration (s)');
    grid(True); xlim((0,12)); ylim((0,5));

#####################################################################################

#####################################################################################
#              slide 19 -- spacing between hotspots

# params for slide 19
def hsSpaceParams (space = [50, 25, 12.5, 6.25]):
  lsec,lopt,lval = [],[],[]
  for b in space:
    lsec.append(["run","run","run","set","set"])
    lopt.append(["simstr","runit","saveout","ip3_origin","boost_every"])
    lval.append(["hsSpaceTest_"+str(b),"1","1",str(2.0*12040.0),str(b)])
  return lsec, lopt, lval

# run sims for slide 19 and save output ( baseline, double IP3R, quadruple ER )
def hsSpaceRun (skip=[]):
  shortRun(hsSpaceParams,skip)

# read results from slide 19 sims
def hsSpaceRead ():
  return shortRead(hsSpaceParams)

# draw results from slide 19
def hsSpaceDraw (dat=None):
  if dat is None: dat = hsSpaceRead()
  from cawave import wavenq,mydraw,data,recdt
  tx, ty = -0.15, 1.06
  ltitles = ["50", "25", "12.5", "6.25"];
  tlx = ["a","b","c","d"];
  for i in xrange(len(ltitles)):
    nq = wavenq(dat[i]["cytca"]); nq.select("overlap",0,"up",1,"y",">=",500); ld=nqs2pyd(nq); nqsdel(nq);
    ax=subplot(2,len(ltitles),i+1); 
    im=imshow(dat[i]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')    
    if i == 0: ylabel('Position (um)');
    ylim((200,800)); xlim((0,25));
    # if gn == len(ltitles): colorbar(im, use_gridspec=True)
    # title("Boost "+ltitles[i]+ " uM : " + r'$Ca^{2+}_{cyt}$ (mM)'); 
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[i],fontsize=14,fontweight='bold',transform=ax.transAxes)
    subplot(2,len(ltitles),i+5);
    plot(ld['startt']/1e3,ld['durt']/1e3,'k');
    plot(ld['startt']/1e3,ld['durt']/1e3,'ko');
    xlabel('Time(s)'); grid(True);
    xlim((0,25)); ylim((0,5));
    if i == 0: ylabel('Duration (s)'); 

#####################################################################################

#####################################################################################
#              -- spacing between coldspots

# params for coldspot spacing
def csSpaceParams (space = [50, 25, 12.5, 6.25]):
  lsec,lopt,lval = [],[],[]
  for b in space:
    lsec.append(["run","run","run","set","set"])
    lopt.append(["simstr","runit","saveout","ip3_origin","boost_every"])
    lval.append(["csSpaceTest_"+str(b),"1","1",str(0.85*12040.0),str(b)])
  return lsec, lopt, lval

# run sims for coldspot spacing and save output
def csSpaceRun (skip=[]):
  shortRun(csSpaceParams,skip)

# read results from coldspot spacing sims
def csSpaceRead ():
  return shortRead(csSpaceParams)

# draw results from coldspot spacing
def csSpaceDraw (dat=None):
  if dat is None: dat = csSpaceRead()
  from cawave import wavenq,mydraw,data,recdt
  tx, ty = -0.15, 1.06
  ltitles = ["50", "25", "12.5", "6.25"];
  tlx = ["a","b","c","d"];
  for i in xrange(len(ltitles)):
    nq = wavenq(dat[i]["cytca"]); nq.select("overlap",0,"up",1,"y",">=",500); ld=nqs2pyd(nq); nqsdel(nq);
    ax=subplot(2,len(ltitles),i+1); 
    im=imshow(dat[i]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')    
    if i == 0: ylabel('Position (um)');
    ylim((400,600)); xlim((0,10));
    # if gn == len(ltitles): colorbar(im, use_gridspec=True)
    # title("Boost "+ltitles[i]+ " uM : " + r'$Ca^{2+}_{cyt}$ (mM)'); 
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[i],fontsize=14,fontweight='bold',transform=ax.transAxes)
    subplot(2,len(ltitles),i+5);
    plot(ld['startt']/1e3,ld['durt']/1e3,'k');
    plot(ld['startt']/1e3,ld['durt']/1e3,'ko');
    xlabel('Time(s)'); grid(True);
    xlim((0,10)); ylim((0,5));
    if i == 0: ylabel('Duration (s)'); 

#####################################################################################


#####################################################################################
#              slide 20 -- hot spots that decrease spread

# params for slide 20
def hsDecParams (b=100,fctr=2.0):
  lsec,lopt,lval = [],[],[]
  # baseline
  lsec.append(["run","run","run"])
  lopt.append(["simstr","runit","saveout"])
  lval.append(["hsDecBase_","1","1"])
  # hotspots which are too closely spaced - causing decreased duration wave
  lsec.append(["run","run","run","set","set"])
  lopt.append(["simstr","runit","saveout","ip3_origin","boost_every"])
  simstr = "hsDecTest_boost_every_"+str(b)+"_ip3_origin_"+str(fctr*12040.0)
  lval.append([simstr,"1","1",str(fctr*12040.0),str(b)])
  return lsec, lopt, lval

# run sims for slide 20 and save output ( baseline, decreased wave spread via hotspot )
def hsDecRun (skip=[]):
  shortRun(hsDecParams,skip)

# read results from slide 20 sims
def hsDecRead ():
  return shortRead(hsDecParams)

# draw results from slide 20
def hsDecDraw (dat=None):
  if dat is None: dat = hsDecRead()
  gn = 1; ltitles = ["Baseline", "Hotspots"];
  for i in xrange(len(ltitles)):
    subplot(len(ltitles),1,gn); 
    if gn == len(ltitles): xlabel('Time(s)');
    imshow(dat[i]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')    
    ylabel('Position (um)');
    colorbar()
    title(ltitles[i]); gn += 1;
  show()

#####################################################################################

#####################################################################################
#              slide 21 -- varies gip3r everywhere

# params for slide 21
def gIP3RParams (lfctr=[0.9, 1.0, 1.1]):
  lsec,lopt,lval = [],[],[]
  for fctr in lfctr:
    lsec.append(["run","run","run","set","set"])
    lopt.append(["simstr","runit","saveout","ip3_origin","ip3_notorigin"])
    simstr = "gIP3RTest_gip3r_"+str(fctr*12040.0)
    lval.append([simstr,"1","1",str(fctr*12040.0),str(fctr*12040.0)])
  return lsec, lopt, lval

# run sims for slide 21 and save output 
def gIP3RRun (skip=[]):
  shortRun(gIP3RParams,skip)

# read results from slide 21 sims
def gIP3RRead ():
  return shortRead(gIP3RParams)

# draw results from slide 21
def gIP3RDraw (dat=None,lfctr=[0.9, 1.0, 1.1]):
  if dat is None: dat = gIP3RRead()
  gn,i = 1,0
  tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c"]
  for fctr in lfctr:
    ax=subplot(1,len(lfctr),gn); xlabel('Time(s)');
    imshow(dat[i]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')    
    xlim((0,15)); ylim((400,600));
    if gn == 1: ylabel('Position (um)');
    # if gn == len(lfctr): colorbar();
    # title("gIP3R: " + str(fctr) + r"$X$: " 
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[i],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1; i += 1;

#####################################################################################

#####################################################################################
#              slide 22,23 -- varies gserca everywhere

# params for slide 22,23
def gSERCAParams (lfctr=[0.4, 0.5, 0.75]):
  lsec,lopt,lval = [],[],[]
  for fctr in lfctr:
    lsec.append(["run","run","run","set"])
    lopt.append(["simstr","runit","saveout","gserca"])
    simstr = "gSERCATest_"+str(fctr*0.3913)
    lval.append([simstr,"1","1",str(fctr*0.3913)])
  return lsec, lopt, lval

# run sims for slide 22,23 and save output 
def gSERCARun (skip=[]):
  shortRun(gSERCAParams,skip)

# read results from slide 22,23 sims
def gSERCARead ():
  return shortRead(gSERCAParams)

# draw results from slide 22,23
def gSERCADraw (dat=None,lfctr=[0.4, 0.5, 0.75]):
  if dat is None: dat = gSERCARead()
  gn,i = 1,0
  tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"]
  for fctr in lfctr:
    ax=subplot(1,len(lfctr),gn); xlabel('Time(s)');
    imshow(dat[i]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((275,725)); 
    if gn == 1: ylabel('Position (um)');
    #if gn == len(lfctr): colorbar();
    #title("gserca: " + str(fctr) + "X"); 
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[i],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1; i += 1;

#####################################################################################

#####################################################################################
#              slide 24 -- long batch varying gserca and gip3

######################
# gip3 variations

#
def gIP3RBatchLims ():
  return numpy.linspace(0.8,2.0,109)

# params for slide 24 - part that varies gip3
def gIP3RBatchParams (lfctr=gIP3RBatchLims()):
  lsec,lopt,lval = [],[],[]
  for fctr in lfctr:
    lsec.append(["run","run","run","set","set","set","set","set","set"])
    lopt.append(["simstr","runit","saveout","ip3_origin","ip3_notorigin","ip3_stim","ip3_stimT","gleak","gserca"])
    simstr = "gIP3RBatch_"+str(fctr*120400.0)
    lval.append([simstr,"1","1",str(fctr*120400.0),str(fctr*120400.0),"1.25","2e3","18.06","1.9565"])
  return lsec, lopt, lval

# run sims for gipr part of slide 24 and save output 
def gIP3RBatchRun (skip=[],blog="gIP3RBatch/gIP3RBatch.log",bdir="gIP3RBatch",qsz=10):
  batchRun(gIP3RBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read results from slide 24 sims
def gIP3RBatchRead ():
  return shortRead(gIP3RBatchParams)

#
def subxgtzero (lon, startt):
  plon = []
  for i in xrange(len(lon)):
    if lon[i] > 0:
      plon.append(lon[i]-startt)
    else:
      plon.append(0)
  return plon

# draw figure showing sensitivity to IP3R density
def gIP3RBatchDraw (dat=None,lnq=None,ldw=None,lsp=None,ldur=None,lamp=None,lon=None,xl=(0.75,1.8),stimT=2000, lidx=[11, 18, 86]):
  if dat is None: dat = gIP3RBatchRead()
  if lnq is None: lnq,ldw,lsp,ldur,lamp,lon=getwavenqs(gIP3RBatchParams,dat)
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"]#; lidx = [15, 86]
  G = gridspec.GridSpec(2, 4)
  for idx in [lidx[0], lidx[-1]]:
    ax=subplot(G[0,gn*2:gn*2+2]);
    xlabel('Time(s)');
    imshow(dat[idx]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((200,800)); xlim((2,7));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
#  lidx = [15, 18, 86] # reset lidx to include baseline
  sz1,sz2 = 10,10; gfctr = gIP3RBatchLims();
  ax=subplot(G[1,0]); xlim(xl); grid(True);
  text(tx,ty,"c",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plon = subxgtzero(lon,stimT)
  plot(gfctr[0:86],plon[0:86],'ko',markersize=sz1); xlabel("IP3R density"); title("Onset (ms)");
  gtt = [gfctr[lidx[0]], gfctr[lidx[1]] ,gfctr[lidx[2]]]; 
  ltt = [plon[lidx[0]], plon[lidx[1]], plon[lidx[2]] ]; plot(gtt,ltt,'ro',markersize=sz2);
  # print the values of red dots
  print '\nred dots values for onset:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nip3rDensity={1}, onset(ms)={2}'.format(val[0],val[1],ltt[val[0]])
  ylim((0,250));
  ax=subplot(G[1,1]); xlim(xl); grid(True);
  text(tx,ty,"d",fontsize=14,fontweight='bold',transform=ax.transAxes);
  #plot(gfctr,lsp,'ko',markersize=sz1); xlabel("IP3R density"); title(r"Speed ($\mu$m/s)");
  plot(gfctr[0:86],lsp[0:86],'ko',markersize=sz1); xlabel("IP3R density"); title(r"Speed ($\mu$m/s)");
  ltt = [lsp[lidx[0]], lsp[lidx[1]], lsp[lidx[2]] ]; plot(gtt,ltt,'ro',markersize=sz2);
  # print the values of red dots
  print '\nred dots values for speed:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nip3rDensity={1}, speed(um/s)={2}'.format(val[0],val[1],ltt[val[0]])
  ax=subplot(G[1,2]); xlim(xl); grid(True);
  text(tx,ty,"e",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr[0:86],numpy.array(ldur)[0:86]/1e3,'ko',markersize=sz1); xlabel("IP3R density"); title("Duration (s)");
  ltt = [ldur[lidx[0]]/1e3, ldur[lidx[1]]/1e3, ldur[lidx[2]]/1e3 ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0,1.2));
  # print the values of red dots
  print '\nred dots values for duration:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nip3rDensity={1}, duration(s)={2}'.format(val[0],val[1],ltt[val[0]])
  ax=subplot(G[1,3]); xlim(xl); grid(True); 
  text(tx,ty,"f",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr[0:86],lamp[0:86],'ko',markersize=sz1); xlabel("IP3R density"); title("Amplitude (mM)");
  ltt = [lamp[lidx[0]], lamp[lidx[1]], lamp[lidx[2]] ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0,0.0021));
  print '\nred dots values for amplitude:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nip3rDensity={1}, amplitude(mM)={2}'.format(val[0],val[1],ltt[val[0]])
  tight_layout()

##############################################################################################################
# test ip3diff (and cadiff??) in the continuous model

#
def IP3DIFFBatchLims ():
  return numpy.linspace(0,10,101)

#
def IP3DIFFBatchParams (ldiff=IP3DIFFBatchLims()):
  lsec,lopt,lval = [],[],[]
  alls,allo,allv = [],[],[] # params shared across sims in this batch
  NewParam(alls,allo,allv,"set","ip3_stim","1.25")
  NewParam(alls,allo,allv,"set","ip3_stimT","2e3")
  NewParam(alls,allo,allv,"set","boost_halfw","5.0")
  NewParam(alls,allo,allv,"set","ip3_origin","120400.0")
  NewParam(alls,allo,allv,"set","ip3_notorigin","120400.0")
  NewParam(alls,allo,allv,"set","gleak",str(18.06))
  NewParam(alls,allo,allv,"set","gserca",str(1.9565))
  NewParam(alls,allo,allv,"run","runit","1")
  NewParam(alls,allo,allv,"run","saveout","1")
  for IP3DIFF in 1.415 * ldiff: # variations in ip3 diffusion coefficient
    tmpsec,tmpopt,tmpval=[],[],[]
    for i in xrange(len(alls)): tmpsec.append(alls[i]); tmpopt.append(allo[i]); tmpval.append(allv[i]);
    simstr = "IP3DIFF_"+str(IP3DIFF)
    NewParam(tmpsec,tmpopt,tmpval,"set","ip3Diff",str(IP3DIFF))
    NewParam(tmpsec,tmpopt,tmpval,"run","simstr",simstr)
    lsec.append(tmpsec); lopt.append(tmpopt); lval.append(tmpval);
  return lsec, lopt, lval

#
def IP3DIFFBatchRun (skip=[],blog="IP3DIFFBatch/IP3DIFFBatch.log",bdir="IP3DIFFBatch",qsz=12):
  batchRun(IP3DIFFBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

#
def IP3DIFFBatchRead ():
  return shortRead(IP3DIFFBatchParams)

#
def IP3DIFFBatchDraw (dat=None,lnq=None,ldw=None,lsp=None,ldur=None,lamp=None,lon=None,xl=(0.12,2.0),stimT=2000, lidx=[1, 10, 14]):
  if dat is None: dat = gIP3RBatchRead()
  if lnq is None: lnq,ldw,lsp,ldur,lamp,lon=getwavenqs(IP3DIFFBatchParams,dat)
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"]#; lidx = [15, 86]
  G = gridspec.GridSpec(2, 4)
  for idx in [lidx[0], lidx[-1]]:
    ax=subplot(G[0,gn*2:gn*2+2]);
    xlabel('Time(s)');
    imshow(dat[idx]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((300,700)); xlim((2,5));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
  sz1,sz2 = 10,10; gfctr = IP3DIFFBatchLims() * 1.415;
  xltxt=r'$IP_3$ diff ($\mu$m$^2$/ms)'
  sidx,eidx=1,15
  ax=subplot(G[1,0]); xlim(xl); grid(True);
  text(tx,ty,"c",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plon = subxgtzero(lon,stimT)
  plot(gfctr[sidx:eidx],plon[sidx:eidx],'ko',markersize=sz1); xlabel(xltxt); title("Onset (ms)");
  gtt = [gfctr[lidx[0]], gfctr[lidx[1]] ,gfctr[lidx[2]]]; 
  ltt = [plon[lidx[0]], plon[lidx[1]], plon[lidx[2]] ]; plot(gtt,ltt,'ro',markersize=sz2);
  # print the values of red dots
  print '\nred dots values for onset:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nIP3DIFF={1}, onset(ms)={2}'.format(val[0],val[1],ltt[val[0]])
  ylim((0,250));
  ax=subplot(G[1,1]); xlim(xl); grid(True);
  text(tx,ty,"d",fontsize=14,fontweight='bold',transform=ax.transAxes);
  #plot(gfctr,lsp,'ko',markersize=sz1); xlabel(xltxt); title(r"Speed ($\mu$m/s)");
  plot(gfctr[sidx:eidx],lsp[sidx:eidx],'ko',markersize=sz1); xlabel(xltxt); title(r"Speed ($\mu$m/s)");
  ltt = [lsp[lidx[0]], lsp[lidx[1]], lsp[lidx[2]] ]; plot(gtt,ltt,'ro',markersize=sz2);
  # print the values of red dots
  print '\nred dots values for speed:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nIP3DIFF={1}, speed(um/s)={2}'.format(val[0],val[1],ltt[val[0]])
  ax=subplot(G[1,2]); xlim(xl); grid(True);
  text(tx,ty,"e",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr[sidx:eidx],numpy.array(ldur)[sidx:eidx]/1e3,'ko',markersize=sz1); xlabel(xltxt); title("Duration (s)");
  ltt = [ldur[lidx[0]]/1e3, ldur[lidx[1]]/1e3, ldur[lidx[2]]/1e3 ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0,1.2));
  # print the values of red dots
  print '\nred dots values for duration:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nIP3DIFF={1}, duration(s)={2}'.format(val[0],val[1],ltt[val[0]])
  ax=subplot(G[1,3]); xlim(xl); grid(True); 
  text(tx,ty,"f",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr[sidx:eidx],lamp[sidx:eidx],'ko',markersize=sz1); xlabel(xltxt); title("Amplitude (mM)");
  ltt = [lamp[lidx[0]], lamp[lidx[1]], lamp[lidx[2]] ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0,0.0021));
  print '\nred dots values for amplitude:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nIP3DIFF={1}, amplitude(mM)={2}'.format(val[0],val[1],ltt[val[0]])
  tight_layout()

##############################################################################################################


######################
# test IP3R hot-spots (modulates IP3R density at hotspots and spacing of hotspots) - HSIP3RBatch

# limits for HS densities
def HSIP3RBatchLimsDense ():
  return numpy.linspace(0.8,2,55)

# limits for HS spacing
def HSIP3RBatchLimsSpace ():
  return numpy.linspace(15,100,18)

# params for batch varying IP3R hotspots (density and spacing)
def HSIP3RBatchParams (lfdense=HSIP3RBatchLimsDense(),lfspace=HSIP3RBatchLimsSpace(),cs=0.8):
  lsec,lopt,lval = [],[],[]
  alls,allo,allv = [],[],[] # params shared across sims in this batch
  NewParam(alls,allo,allv,"set","ip3_stim","1.25")
  NewParam(alls,allo,allv,"set","ip3_stimT","2e3")
  NewParam(alls,allo,allv,"set","gleak","18.06")
  NewParam(alls,allo,allv,"set","gserca","1.9565")
  NewParam(alls,allo,allv,"set","boost_halfw","5.0")
  NewParam(alls,allo,allv,"set","ip3_notorigin",str(cs*120400.0))
  NewParam(alls,allo,allv,"run","tstop","15000")
  for dense in lfdense: # variations in density of IP3R hotspots
    for space in lfspace: # variations in distance between IP3R hotspots
      tmpsec,tmpopt,tmpval=[],[],[]
      for i in xrange(len(alls)): tmpsec.append(alls[i]); tmpopt.append(allo[i]); tmpval.append(allv[i]);
      simstr = "HSIP3RBatch_Dense_"+str(dense*120400.0)+"_Space_"+str(space)
      NewParam(tmpsec,tmpopt,tmpval,"run","simstr",simstr)
      NewParam(tmpsec,tmpopt,tmpval,"run","runit","1")
      NewParam(tmpsec,tmpopt,tmpval,"run","saveout","1")
      NewParam(tmpsec,tmpopt,tmpval,"set","ip3_origin",str(dense*120400.0))
      NewParam(tmpsec,tmpopt,tmpval,"set","boost_every",str(space))
      lsec.append(tmpsec); lopt.append(tmpopt); lval.append(tmpval);
  return lsec, lopt, lval

# run this batch
def HSIP3RBatchRun (skip=[],blog="HSIP3RBatch/HSIP3RBatch.log",bdir="HSIP3RBatch",qsz=12):
  batchRun(HSIP3RBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read data from this batch
def HSIP3RBatchRead ():
  return shortRead(HSIP3RBatchParams)

# print the values at the row in the nqs (div by baseline IP3R density of 120400.0)
def HSIP3RBatchParamsAtRow (nq,row):
  nq.tog("DB")
  simstr = nq.get("simstr",row).s
  dense = nq.getcol("ip3_origin").x[row] / 120400.0
  cs = nq.getcol("ip3_notorigin").x[row] / 120400.0
  space = nq.getcol("boost_every").x[row] 
  print "row " , row, ": dense: " , dense, " space : " , space, " CS: " , cs, ", simstr:",simstr
  speed, amp = nq.getcol("speed").x[row],nq.getcol("amp").x[row]  
  dur, onset = nq.getcol("dur").x[row],nq.getcol("onset").x[row]  
  print "\t speed:",speed,", amp:",amp,", dur:",dur,", onset:",onset
  return dense,cs,space,speed,amp,dur,onset

# draw the waves for first part of the figure
def HSIP3RBatchDrawWaves (nq,lidx=[163, 577, 973, 865, 871, 877]):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"];
  lfdense,lfspace = HSIP3RBatchLimsDense(),HSIP3RBatchLimsSpace()
  nrows,ncols = 2,3
  for idx in lidx:
    ax=subplot(nrows,ncols,gn+1);
    dat = loadfromnq(nq,idx);
    HSIP3RBatchParamsAtRow(nq,idx) # print info
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((0,1000)); xlim((2,10));
    if gn == 0 or gn == 3: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    if gn > 2: xlabel('Time (s)');
    gn += 1;  

# draw waves/features 
def HSIP3RBatchDrawSub (nq,G,gRow,fixedspace=True):
  lfdense,lfspace = HSIP3RBatchLimsDense(),HSIP3RBatchLimsSpace()
  lidx,gctr,lidx2,xtxt=[],[],[],"";
  if fixedspace: # fixed spacing
    lidx = [109, 973] # index into nqs # new run of sims
    #lidx = [163, 973] # index into nqs
    nq.select("boost_every",20.0)
    xtxt = r"Density"
    gfctr = lfdense
  else: # fixed density
    lidx = [864, 881] # index into nqs # new run of sims
    #lidx = [865, 877] # index into nqs
    nq.select("ip3_origin",lfdense[48]*120400.0)
    xtxt = r"Spacing ($\mu$m)"
    gfctr = lfspace
  lsp,ldur,lamp,lon,ltmp = getwpropcols(nq,db=False); nq.tog("DB")
  #return ltmp # dubgging code
  lidx2 = [ltmp.index(lidx[0]), ltmp.index(lidx[1])]; 
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c"]
  if not fixedspace: tlx = ["d", "e", "f"]; 
  for idx in lidx: # draw the waves
    dat = loadfromnq(nq,idx); HSIP3RBatchParamsAtRow(nq,idx) # load data, print info
    ax=subplot(G[gRow,gn:gn+1]); xlabel('Time(s)');
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,15,0,1000),origin='lower')
    ylim((300,700)); xlim((2,6));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
  ly = [ (0,115) ]; ltitles = [r"Speed ($\mu$m/s)"]; gn = 2; sz1,sz2=10,10;  
  scaley = [1.0,1.0,1.0/1e3,1.0];  naxbins=6; lvals=[ lsp ];
  for pdx in xrange(1): # plots for the four wave features
    ax=subplot(G[gRow,gn:gn+1]); xlabel(xtxt); title(ltitles[pdx]); ylim(ly[pdx]);
    xlim([amin(gfctr)-amin(gfctr)/5.0,amax(gfctr)+amin(gfctr)/5.0]); # to avoid data points being flush with graph borders
    #xlim([amin(gfctr)-amin(gfctr)/10.0,amax(gfctr)+amin(gfctr)/10.0]);
    ax.locator_params(nbins=naxbins);
    text(tx,ty,tlx[pdx+2],fontsize=14,fontweight='bold',transform=ax.transAxes); grid(True);
    plot(gfctr,numpy.array(lvals[pdx])*scaley[pdx],'ko',markersize=sz1); 
    gtt,ltt=[gfctr[lidx2[0]],gfctr[lidx2[1]]],[lvals[pdx][lidx2[0]],lvals[pdx][lidx2[1]]];
    plot(gtt,numpy.array(ltt)*scaley[pdx],'ro',markersize=sz1);
    print '\nred dots values for ' + xtxt + ' :'
    for val in enumerate(gtt):
      print 'red dot index({0}):\n{3}={1}, Speed(um)={2}'.format(val[0],val[1],ltt[val[0]],xtxt)
    gn += 1;


# draw waves and speeds
def HSIP3RBatchDraw (nq):
  nrows,ncols = 2,3
  G = gridspec.GridSpec(nrows, ncols);
  HSIP3RBatchDrawSub(nq,G,0,fixedspace=True) # fixed spacing
  HSIP3RBatchDrawSub(nq,G,1,fixedspace=False) # fixed density
  tight_layout()


# draw the wave properties stored in the nqs (from HSIP3RBatchNQS)
def HSIP3RBatchDrawWaveProps (nq,xl=(15,100),yl=(0.8,2),stimt=2e3):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"]; 
  lcols = ["speed", "onset"]
  ltitles = [r'speed ($\mu$m/s)','onset (ms)'];
  lvmin,lvmax = [40, 0, 0], [80, 300, 0.002]
  for col in lcols:
    ax=subplot(1,len(lcols),gn+1); S,extent = getmatbyparams(nq,col);
    if col == "time":
      imshow(S/1e3,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=extent,origin='lower',interpolation='None')
    elif col == "onset":
      imshow(S-stimt,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=extent,origin='lower',interpolation='None')
    else:
      imshow(S,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=extent,origin='lower',interpolation='None')
    ylim(yl); xlim(xl); xlabel(r'Hotspot spacing ($\mu$m)'); title(ltitles[gn]);
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes); colorbar();
    if gn == 0: ylabel(r'IP3R density');    
    gn += 1

# get the NQS with batch data / wave properties
def HSIP3RBatchNQS ():
  return NQS("data/13sep09_HSIP3RBatch_wp2.nqs")

# creates a matrix with values arranged by HS,CS values
# nq from params2nq -> addwavenqcol -> addwavepropcols
# col specifies which column from nqs to use
def getmatbyparams (nq,colval,xparm=HSIP3RBatchLimsSpace,yparm=HSIP3RBatchLimsDense,coly="ip3_origin",colx="boost_every",scaley=120400.0,scalex=1):
  lx,ly = xparm(),yparm()
  S = numpy.zeros((len(ly),len(lx)))
  for i in xrange(len(ly)):
    for j in xrange(len(lx)):
      N = nq.select(-1,coly,scaley*ly[i],colx,scalex*lx[j])
      if N < 1: continue
      idx = int(nq.ind.x[0])
      S[i][j] = nq.getcol(colval).x[idx]
  extent=amin(lx),amax(lx),amin(ly),amax(ly)
  return S,extent

######################
# test ER stacks (density and spacing) - STACKERBatch

# limits for STACK densities
def STACKERBatchLimsDense ():
  return numpy.linspace(0.8,2.0,55)

# limits for STACK spacing
def STACKERBatchLimsSpace ():
  return numpy.linspace(15,100,18)

# get the NQS with batch data / wave properties
def STACKERBatchNQS ():
  return NQS("data/13sep11_STACKERBatch_wp2.nqs")

# params for batch varying ER STACKS 
#  cs is density between the stacks
def STACKERBatchParams (lfdense=STACKERBatchLimsDense(),lfspace=STACKERBatchLimsSpace(),cs=0.8):
  lsec,lopt,lval = [],[],[]
  alls,allo,allv = [],[],[] # params shared across sims in this batch
  NewParam(alls,allo,allv,"set","ip3_stim","1.25")
  NewParam(alls,allo,allv,"set","ip3_stimT","2e3")
  NewParam(alls,allo,allv,"set","boost_halfw","5.0")
  gip3r, gsercar, gleakr = cs*120400.0, cs*1.9565, cs*18.06
  NewParam(alls,allo,allv,"set","ip3_origin",str(gip3r))
  NewParam(alls,allo,allv,"set","ip3_notorigin",str(gip3r))
  NewParam(alls,allo,allv,"set","gleak",str(gleakr))
  NewParam(alls,allo,allv,"set","gserca",str(gsercar))
  NewParam(alls,allo,allv,"run","runit","1")
  NewParam(alls,allo,allv,"run","saveout","1")
  allsimstr = "STACKERBatch_BoostHW5_IP3R_" + str(gip3r) + "_GLEAK_"+str(gleakr)+"_GSERCA_"+str(gsercar)
  for dense in lfdense: # variations in density of ER stacks
    for space in lfspace: # variations in distance between ER stacks
      tmpsec,tmpopt,tmpval=[],[],[]
      for i in xrange(len(alls)): tmpsec.append(alls[i]); tmpopt.append(allo[i]); tmpval.append(allv[i]);
      ers = dense*1.0/cs
      simstr = allsimstr + "_DENSE_" + str(dense) + "_ERScale_"+str(ers)+"_Space_"+str(space)
      NewParam(tmpsec,tmpopt,tmpval,"run","simstr",simstr)
      NewParam(tmpsec,tmpopt,tmpval,"set","er_scale",str(ers))
      NewParam(tmpsec,tmpopt,tmpval,"set","boost_every",str(space))
      lsec.append(tmpsec); lopt.append(tmpopt); lval.append(tmpval);
  return lsec, lopt, lval

# run this batch
def STACKERBatchRun (skip=[],blog="STACKERBatch/STACKERBatch.log",bdir="STACKERBatch",qsz=12):
  batchRun(STACKERBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read data from this batch
def STACKERBatchRead ():
  return shortRead(STACKERBatchParams)

# print the values at the row in the nqs (div by baseline ER density)
def STACKERBatchParamsAtRow (nq,row):
  nq.tog("DB")
  simstr = nq.get("simstr",row).s
  gip3r = nq.getcol("ip3_origin").x[row] / 120400.0
  cs = nq.getcol("ip3_notorigin").x[row] / 120400.0
  space = nq.getcol("boost_every").x[row] 
  ers = nq.getcol("er_scale").x[row]
  gleak = nq.getcol("gleak").x[row]
  gserca = nq.getcol("gserca").x[row]
  print "row ",row,",gip3r:",gip3r,",space:",space,",CS:",cs,"er_scale:",ers,",gleak:",gleak,",gserca:",gserca,",simstr:",simstr
  speed, amp = nq.getcol("speed").x[row],nq.getcol("amp").x[row]  
  dur, onset = nq.getcol("dur").x[row],nq.getcol("onset").x[row]  
  print "\t speed:",speed,", amp:",amp,", dur:",dur,", onset:",onset
  return gip3r,cs,space,ers,gleak,gserca,row,simstr

# draw the waves for first part of the figure
def STACKERBatchDrawWaves (nq,lidx=[163, 577, 973, 865, 871, 877]):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"];
  lfdense,lfspace = STACKERBatchLimsDense(),STACKERBatchLimsSpace()
  nrows,ncols = 2,3
  for idx in lidx:
    ax=subplot(nrows,ncols,gn+1);
    dat = loadfromnq(nq,idx);
    STACKERBatchParamsAtRow(nq,idx) # print info
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((0,1000)); xlim((2,10));
    if gn == 0 or gn == 3: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    if gn > 2: xlabel('Time (s)');
    gn += 1;  

# draw the wave properties stored in the nqs (from STACKERBatchNQS)
def STACKERBatchDrawWaveProps  (nq,xl=(15,100),yl=(0.8,2),stimt=2e3,cs=0.8):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"]; 
  lcols = ["speed", "onset"]
  ltitles = [r'speed ($\mu$m/s)','onset (ms)'];
  lvmin,lvmax = [40, 0, 0], [80, 300, 0.002]
  for col in lcols:
    ax=subplot(1,len(lcols),gn+1);
    S,extent = getmatbyparams(nq,col,xparm=STACKERBatchLimsSpace,yparm=STACKERBatchLimsDense,coly="er_scale",colx="boost_every",scaley=1.0/cs,scalex=1.0);
    if col == "time":
      imshow(S/1e3,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=extent,origin='lower',interpolation='None')
    elif col == "onset":
      imshow(S-stimt,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=extent,origin='lower',interpolation='None')
    else:
      imshow(S,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=extent,origin='lower',interpolation='None')
    ylim(yl); xlim(xl); xlabel(r'Stack spacing ($\mu$m)'); title(ltitles[gn]);
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes); colorbar();
    if gn == 0: ylabel(r'ER density');    
    gn += 1

# draw waves/features 
def STACKERBatchDrawSub (nq,G,gRow,fixedspace=True):
  lfdense,lfspace = STACKERBatchLimsDense(),STACKERBatchLimsSpace()
  lidx,gctr,lidx2,xtxt=[],[],[],"";
  if fixedspace: # fixed spacing
    lidx = [1, 973] # index into nqs
    nq.select("boost_every",20.0)
    xtxt = r"Density"
    gfctr = lfdense
  else: # fixed density
    lidx = [864, 881] # index into nqs 
    nq.select("er_scale",lfdense[48]*1.0/0.8)
    xtxt = r"Spacing ($\mu$m)"
    gfctr = lfspace
  lsp,ldur,lamp,lon,ltmp = getwpropcols(nq,db=False); nq.tog("DB")
  lidx2 = [ltmp.index(lidx[0]), ltmp.index(lidx[1])]; 
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c"]
  if not fixedspace: tlx = ["d", "e", "f"]; 
  for idx in lidx: # draw the waves
    dat = loadfromnq(nq,idx); STACKERBatchParamsAtRow(nq,idx) # load data, print info
    ax=subplot(G[gRow,gn:gn+1]); xlabel('Time(s)');
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((300,700)); xlim((2,6));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
  #ly = [ (65,95) ]; ltitles = [r"Speed ($\mu$m/s)"]; gn = 2; sz1,sz2=10,10;
  ly = [ (0,115) ]; ltitles = [r"Speed ($\mu$m/s)"]; gn = 2; sz1,sz2=10,10; #changed ylimits to be similar to the hotspots graph
  scaley = [1.0,1.0,1.0/1e3,1.0];  naxbins=6; lvals=[ lsp ];
  for pdx in xrange(1): # for the four wave features
    ax=subplot(G[gRow,gn:gn+1]); xlabel(xtxt); title(ltitles[pdx]); ylim(ly[pdx]);
    xlim([amin(gfctr)-amin(gfctr)/10.0,amax(gfctr)+amin(gfctr)/10.0]);
    ax.locator_params(nbins=naxbins);
    text(tx,ty,tlx[pdx+2],fontsize=14,fontweight='bold',transform=ax.transAxes); grid(True);
    plot(gfctr,numpy.array(lvals[pdx])*scaley[pdx],'ko',markersize=sz1); 
    gtt,ltt=[gfctr[lidx2[0]],gfctr[lidx2[1]]],[lvals[pdx][lidx2[0]],lvals[pdx][lidx2[1]]];
    plot(gtt,numpy.array(ltt)*scaley[pdx],'ro',markersize=sz1); 
    gn += 1;

# draw waves and speeds
def STACKERBatchDraw (nq):
  nrows,ncols = 2,3
  G = gridspec.GridSpec(nrows, ncols);
  STACKERBatchDrawSub(nq,G,0,fixedspace=True) # fixed spacing
  STACKERBatchDrawSub(nq,G,1,fixedspace=False) # fixed density
  tight_layout()

######################
# test diffusion coefficient effect on IP3R hotspots

#
def HSDIFFBatchLims ():
  return numpy.linspace(0,10,101)

# params for batch varying Ca2+ diffusion coefficient with the IP3R hotspots
#  hs is IP3R density at hotspots, cs is IP3R density between the hotspots 
def HSDIFFBatchParams (ldiff=HSDIFFBatchLims(),hs=1.75,cs=0.8):
  lsec,lopt,lval = [],[],[]
  alls,allo,allv = [],[],[] # params shared across sims in this batch
  NewParam(alls,allo,allv,"set","ip3_stim","1.25")
  NewParam(alls,allo,allv,"set","ip3_stimT","2e3")
  NewParam(alls,allo,allv,"set","boost_halfw","5.0")
  hsip3r,csip3r = 120400.0*hs,120400.0*cs
  NewParam(alls,allo,allv,"set","ip3_origin",str(hsip3r))
  NewParam(alls,allo,allv,"set","ip3_notorigin",str(csip3r))
  NewParam(alls,allo,allv,"set","gleak",str(18.06))
  NewParam(alls,allo,allv,"set","gserca",str(1.9565))
  NewParam(alls,allo,allv,"run","runit","1")
  NewParam(alls,allo,allv,"run","saveout","1")
  NewParam(alls,allo,allv,"set","boost_every","20.0")
  NewParam(alls,allo,allv,"set","ip3Diff","1.415")
  NewParam(alls,allo,allv,"run","tstop","15000")
  allsimstr = "HSDIFFBatch_BoostEvery20_BoostHW5_HSIP3R_"+str(hsip3r)+"_CSIP3R_"+str(csip3r)
  for caDiff in 0.080 * ldiff: # variations in Ca2+ diffusion coefficient
    tmpsec,tmpopt,tmpval=[],[],[]
    for i in xrange(len(alls)): tmpsec.append(alls[i]); tmpopt.append(allo[i]); tmpval.append(allv[i]);
    simstr = allsimstr + "_caDiff_"+str(caDiff)
    NewParam(tmpsec,tmpopt,tmpval,"run","simstr",simstr)
    NewParam(tmpsec,tmpopt,tmpval,"set","caDiff",str(caDiff))
    lsec.append(tmpsec); lopt.append(tmpopt); lval.append(tmpval);
  return lsec, lopt, lval

# run this batch
def HSDIFFBatchRun (skip=[],blog="HSDIFFBatch/HSDIFFBatch.log",bdir="HSDIFFBatch",qsz=12):
  batchRun(HSDIFFBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read data from this batch
def HSDIFFBatchRead ():
  return shortRead(HSDIFFBatchParams)

# get the NQS with batch data / wave properties
def HSDIFFBatchNQS ():
  return NQS("data/14sep18_HSDIFFBatch_wp2.nqs")

# print the values at the row in the nqs (div by baseline ER density)
def HSDIFFBatchParamsAtRow (nq,row):
  nq.tog("DB")
  simstr = nq.get("simstr",row).s
  ip3Diff = nq.getcol("ip3Diff").x[row]
  caDiff = nq.getcol("caDiff").x[row]
  print "row ",row,"ip3Diff:",ip3Diff,",caDiff:",caDiff,"simstr:",simstr
  speed, amp = nq.getcol("speed").x[row],nq.getcol("amp").x[row]  
  dur, onset = nq.getcol("dur").x[row],nq.getcol("onset").x[row]  
  print "\t speed:",speed,", amp:",amp,", dur:",dur,", onset:",onset
  return ip3Diff,caDiff,row,simstr

# draw the waves for first part of the figure
def HSDIFFBatchDrawWaves (nq,lidx=[ 25*26+0, 25*26+26/2, 25*26+25, 0*26+0,  0*26+26/2,  0*26+25]):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"];
  ldiff = HSDIFFBatchLims()
  nrows,ncols = 2,3
  for idx in lidx:
    ax=subplot(nrows,ncols,gn+1);
    dat = loadfromnq(nq,idx);
    HSDIFFBatchParamsAtRow(nq,idx) # print info
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((500-300,500+300)); xlim((2,6));
    if gn == 0 or gn == 3: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    if tlx[gn] in ["a", "b", "c"]: text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    if gn > 2: xlabel('Time (s)');
    gn += 1;  

# draw the wave properties stored in the nqs (from HSDIFFBatchNQS)
def HSDIFFBatchDrawWaveProps  (nq,xl=(0,5),yl=(0,5),stimt=2e3,cs=0.8):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"]; 
  lcols = ["speed", "onset"]
  ltitles = [r'speed ($\mu$m / s)','onset (ms)',"dist","amp","dur"];
  lvmin,lvmax = [40, 0, 500, 0, 0], [160, 205, 1000, 0.002, 30]
  for col in lcols:
    ax=subplot(1,len(lcols),gn+1);
    S,extent = getmatbyparams(nq,col,xparm=HSDIFFBatchLims,yparm=HSDIFFBatchLims,coly="caDiff",colx="ip3Diff",scaley=0.08,scalex=1.415);
    ext = (0, 1.415*5, 0, 0.08*5)
    if col == "time":
      imshow(S/1e3,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=ext,origin='lower',interpolation='None')
    elif col == "onset":
      imshow(S-stimt,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=ext,origin='lower',interpolation='None')
    else:
      imshow(S,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=ext,origin='lower',interpolation='None')
    # ylim(yl); xlim(xl); 
    xlabel(r'IP3 diffusion coefficent ($\mu$m$^2$ / ms)'); title(ltitles[gn]);
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes); colorbar();
    if gn == 0: ylabel(r'Ca$^{2+}$ diffusion coefficient ($\mu$m$^2$ / ms)');    
    gn += 1

# get the wave property columns
def getwpropcols (nq,db=True):
  lsp,ldur,lamp,lon,lidx = [],[],[],[],[]
  if db:
    nq.tog("DB")
    lidx = numpy.linspace(0,nq.v[0].size()-1,nq.v[0].size())
  lsp = nq.getcol("speed").to_python()
  ldur = nq.getcol("dur").to_python()
  lamp = nq.getcol("amp").to_python()
  lon = nq.getcol("onset").to_python()
  lidx = nq.ind.to_python()
  return lsp,ldur,lamp,lon,lidx

# draw figure showing response to ca2+ diff coeff changes
def HSDIFFBatchDraw (nq,xl=(0,0.8),stimT=2000,lidx=[1,100]):
  lsp,ldur,lamp,lon, returnedlidx = getwpropcols(nq)
  #lsp,ldur,lamp,lon = getwpropcols(nq)
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"];
  G = gridspec.GridSpec(2, 4)
  for idx in lidx:
    ax=subplot(G[0,gn*2:gn*2+2]);
    xlabel('Time(s)');
    dat = loadfromnq(nq,idx);
    HSDIFFBatchParamsAtRow(nq,idx)
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,15,0,1000),origin='lower')
    ylim((0,1000)); xlim((2,5));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
  naxbins = 5
  xtxt = r'Ca$^{2+}$ diff. ($\mu$m$^2$/ms)'
  # lidx = [20, 100] # reset lidx to include baseline
  sz1,sz2 = 10,10; gfctr = 0.08 * HSDIFFBatchLims();
  ax=subplot(G[1,0]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"c",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plon = subxgtzero(lon,stimT)
  plot(gfctr,plon,'ko',markersize=sz1); xlabel(xtxt); title("Onset (ms)");
  gtt = [gfctr[lidx[0]], gfctr[lidx[1]] ];
  ltt = [plon[lidx[0]], plon[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0,50));
  ax=subplot(G[1,1]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"d",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,lsp,'ko',markersize=sz1); xlabel(xtxt); title(r"Speed ($\mu$m/s)");
  ltt = [lsp[lidx[0]], lsp[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0,300))
  ax=subplot(G[1,2]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"e",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,numpy.array(ldur)/1e3,'ko',markersize=sz1); xlabel(xtxt); title("Duration (s)");
  ltt = [ldur[lidx[0]]/1e3, ldur[lidx[1]]/1e3 ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0,1.2));
  ax=subplot(G[1,3]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"f",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,lamp,'ko',markersize=sz1); xlabel(xtxt); title("Amplitude (mM)");
  ltt = [lamp[lidx[0]], lamp[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0.0014,0.0020));
  tight_layout()

######################
# test diffusion coefficient effect on ER stacks

#
def STACKDIFFBatchLims ():
  return numpy.linspace(0,10,101)

# params for batch varying diffusion coefficients with the ER stacks
#  hs is ER density at stacks, cs is ER density between the stacks
def STACKDIFFBatchParams (ldiff=STACKDIFFBatchLims(),hs=1.75,cs=0.8):
  lsec,lopt,lval = [],[],[]
  alls,allo,allv = [],[],[] # params shared across sims in this batch
  NewParam(alls,allo,allv,"set","ip3_stim","1.25")
  NewParam(alls,allo,allv,"set","ip3_stimT","2e3")
  NewParam(alls,allo,allv,"set","boost_halfw","5.0")
  gip3r, gsercar, gleakr = cs*120400.0, cs*1.9565, cs*18.06
  NewParam(alls,allo,allv,"set","ip3_origin",str(gip3r))
  NewParam(alls,allo,allv,"set","ip3_notorigin",str(gip3r))
  NewParam(alls,allo,allv,"set","gleak",str(gleakr))
  NewParam(alls,allo,allv,"set","gserca",str(gsercar))
  ers = hs*1.0/cs
  NewParam(alls,allo,allv,"set","er_scale",str(ers))
  NewParam(alls,allo,allv,"run","runit","1")
  NewParam(alls,allo,allv,"run","saveout","1")
  NewParam(alls,allo,allv,"set","boost_every","20.0")
  NewParam(alls,allo,allv,"set","ip3Diff","1.415")
  allsimstr = "STACKDIFFBatch_BoostEvery20_BoostHW5_gIP3R_"+str(gip3r)+"_er_scale"+str(ers)+"_hs_"+str(hs)+"_cs_"+str(cs)
  for caDiff in 0.080 * ldiff: # variations in Ca diffusion coefficient
    tmpsec,tmpopt,tmpval=[],[],[]
    for i in xrange(len(alls)): tmpsec.append(alls[i]); tmpopt.append(allo[i]); tmpval.append(allv[i]);
    simstr = allsimstr + "_caDiff_"+str(caDiff)
    NewParam(tmpsec,tmpopt,tmpval,"run","simstr",simstr)
    NewParam(tmpsec,tmpopt,tmpval,"set","caDiff",str(caDiff))
    lsec.append(tmpsec); lopt.append(tmpopt); lval.append(tmpval);
  return lsec, lopt, lval

# run this batch
def STACKDIFFBatchRun (skip=[],blog="STACKDIFFBatch/STACKDIFFBatch.log",bdir="STACKDIFFBatch",qsz=12):
  batchRun(STACKDIFFBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read data from this batch
def STACKDIFFBatchRead ():
  return shortRead(STACKDIFFBatchParams)

# get the NQS with batch data / wave properties
def STACKDIFFBatchNQS ():
  return NQS("data/13sep16_STACKDIFFBatch_wp3.nqs")

# print the values at the row in the nqs 
def STACKDIFFBatchParamsAtRow (nq,row):
  nq.tog("DB")
  simstr = nq.get("simstr",row).s
  ip3Diff = nq.getcol("ip3Diff").x[row]
  caDiff = nq.getcol("caDiff").x[row]
  print "row ",row,"ip3Diff:",ip3Diff,",caDiff:",caDiff,"simstr:",simstr
  # code added from HSDIFFBatchParamsAtRow
  #debug_here()# explore the headers of nq
  speed, amp = nq.getcol("speed").x[row],nq.getcol("amp").x[row]  
  dur, onset = nq.getcol("dur").x[row],nq.getcol("onset").x[row]  
  print "\t speed:",speed,", amp:",amp,", dur:",dur,", onset:",onset
  return ip3Diff,caDiff,row,simstr

# draw the waves for first part of the figure
def STACKDIFFBatchDrawWaves (nq,lidx=[0,1,2,3,4,5]):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"];
  ldiff = STACKDIFFBatchLims()
  nrows,ncols = 2,3
  for idx in lidx:
    ax=subplot(nrows,ncols,gn+1);
    dat = loadfromnq(nq,idx);
    STACKDIFFBatchParamsAtRow(nq,idx) # print info
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((0,1000)); xlim((2,10));
    if gn == 0 or gn == 3: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    if gn > 2: xlabel('Time (s)');
    gn += 1;  

# draw the wave properties stored in the nqs (from STACKDIFFBatchNQS)
def STACKDIFFBatchDrawWaveProps  (nq,xl=(0,5),yl=(0,5),stimt=2e3,cs=0.8):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"]; 
  lcols = ["speed", "onset","dist","amp","dur"]
  ltitles = [r'speed ($\mu$m/s)','onset (ms)',"dist","amp","dur"];
  lvmin,lvmax = [40, 0, 500, 0, 0], [600, 300, 1000, 0.002, 30]
  for col in lcols:
    ax=subplot(1,len(lcols),gn+1);
    S,extent = getmatbyparams(nq,col,xparm=STACKDIFFBatchLims,yparm=STACKDIFFBatchLims,coly="caDiff",colx="ip3Diff",scaley=1.0,scalex=1.0);
    if col == "time":
      imshow(S/1e3,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=extent,origin='lower',interpolation='None')
    elif col == "onset":
      imshow(S-stimt,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=extent,origin='lower',interpolation='None')
    else:
      imshow(S,vmin=lvmin[gn],vmax=lvmax[gn],aspect='auto',extent=extent,origin='lower',interpolation='None')
    ylim(yl); xlim(xl); xlabel(r'IP3 diffusion coefficent'); title(ltitles[gn]);
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes); colorbar();
    if gn == 0: ylabel(r'Ca$^{2+}$ diffusion coefficient');    
    gn += 1

# draw figure showing response to ca2+ diff coeff changes
#def STACKDIFFBatchDraw (nq,xl=(0,0.8),stimT=2000,lidx=[4,100]):
def STACKDIFFBatchDraw (nq,xl=(0,0.8),stimT=2000,lidx=[0,100]):
  lsp,ldur,lamp,lon, returnedlidx = getwpropcols(nq)
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"];
  G = gridspec.GridSpec(2, 4)
  for idx in lidx:
    ax=subplot(G[0,gn*2:gn*2+2]);
    xlabel('Time(s)');
    dat = loadfromnq(nq,idx);
    STACKDIFFBatchParamsAtRow(nq, idx)
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((0,1000)); xlim((2,5));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
  naxbins = 5
  xtxt = r'Ca$^{2+}$ diff. ($\mu$m$^2$/ms)'
  # lidx = [20, 100] # reset lidx to include baseline
  sz1,sz2 = 10,10; gfctr = 0.08 * STACKDIFFBatchLims();
  ax=subplot(G[1,0]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"c",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plon = subxgtzero(lon,stimT)
  plot(gfctr,plon,'ko',markersize=sz1); xlabel(xtxt); title("Onset (ms)");
  gtt = [gfctr[lidx[0]], gfctr[lidx[1]] ];
  ltt = [plon[lidx[0]], plon[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((30,75));
  ax=subplot(G[1,1]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"d",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,lsp,'ko',markersize=sz1); xlabel(xtxt); title(r"Speed ($\mu$m/s)");
  ltt = [lsp[lidx[0]], lsp[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0,300))
  ax=subplot(G[1,2]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"e",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,numpy.array(ldur)/1e3,'ko',markersize=sz1); xlabel(xtxt); title("Duration (s)");
  ltt = [ldur[lidx[0]]/1e3, ldur[lidx[1]]/1e3 ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0.5,1.0));
  ax=subplot(G[1,3]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"f",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,lamp,'ko',markersize=sz1); xlabel(xtxt); title("Amplitude (mM)");
  ltt = [lamp[lidx[0]], lamp[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=sz2); ylim((0.0012,0.0021));
  tight_layout()

######################
# test hot/cold-spots (ER vs IP3R changes) - HSCSIP3RERBatch

#
def HSCSIP3RERBatchLims ():
  return numpy.linspace(0,5,51)

# params for hot/cold-spot variations in IP3R and ER levels
def HSCSIP3RERBatchParams (lfctr=HSCSIP3RERBatchLims()):
  lsec,lopt,lval = [],[],[]
  for fctr in lfctr: # first, the variations in IP3R
    lsec.append(["run","run","run","set","set"])
    lopt.append(["simstr","runit","saveout","ip3_origin","boost_every"])
    simstr = "HSCSIP3RERBatch_Boost50_IP3R_"+str(fctr*120400.0)
    lval.append([simstr,"1","1",str(fctr*120400.0), "50.0"])
  for fctr in lfctr: # second, the variations in ER
    lsec.append(["run","run","run","set","set"])
    lopt.append(["simstr","runit","saveout","er_scale","boost_every"])
    simstr = "HSCSIP3RERBatch_Boost50_ER_"+str(fctr)
    lval.append([simstr,"1","1",str(fctr), "50.0"])
  return lsec, lopt, lval

# run this batch
def HSCSIP3RERBatchRun (skip=[],blog="HSCSIP3RERBatch/HSCSIP3RERBatch.log",bdir="HSCSIP3RERBatch",qsz=12):
  batchRun(HSCSIP3RERBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read data from this batch
def HSCSIP3RERBatchRead ():
  return shortRead(HSCSIP3RERBatchParams)

# draw data from this batch
def HSCSIP3RERBatchDraw (dat=None,lnq=None,ldw=None,lsp=None,ldur=None,lamp=None,lon=None,xl=(0.05,2.25)):
  if dat is None: dat = HSCSIP3RERBatchRead()
  if lnq is None: lnq,ldw,lsp,ldur,lamp,lon=getwavenqs(HSCSIP3RERBatchParams,dat)
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f", "g", "h"]; lidx = [1, 22, 1+51, 22+51];
  gfctr = HSCSIP3RERBatchLims()
  print "a:",lidx[0],gfctr[lidx[0]],"b:",lidx[1],gfctr[lidx[1]]
  G = gridspec.GridSpec(2, 8);
  for idx in lidx:
    ax=subplot(G[0,gn*2:gn*2+2]); xlabel('Time(s)');
    imshow(dat[idx]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((150,850)); xlim((0,8));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
  ly = [ (0,18), (44,57), (1.0,1.06), (0.0018,0.0023) ];
  slc = [ [0,51] , [51, 102] ]; sty1 = ['b2', 'r3']; sty2 = ['bd', 'rd'];  
  ltitles = ["Onset (ms)", r"Speed ($\mu$m/s)", "Duration (s)", "Amplitude (mM)"];
  lletters = ["e", "f", "g", "h"]; lvals = [ lon, lsp, ldur, lamp ]; gn = 0; sz1,sz2=10,10;
  scaley = [1.0,1.0,1.0/1e3,1.0];
  for pdx in xrange(4):
    ax=subplot(G[1,gn:gn+2]); xlim(xl); xlabel("Branch Strength"); title(ltitles[pdx]); ylim(ly[pdx]);
    text(tx,ty,lletters[pdx],fontsize=14,fontweight='bold',transform=ax.transAxes); grid(True);
    for i in xrange(2):
      plot(gfctr,numpy.array(lvals[pdx][slc[i][0]:slc[i][1]])*scaley[pdx],sty1[i],markersize=sz1); 
      gtt,ltt=[gfctr[lidx[0]],gfctr[lidx[1]]],[lvals[pdx][lidx[2*i]],lvals[pdx][lidx[2*i+1]]];
      plot(gtt,numpy.array(ltt)*scaley[pdx],sty2[i],markersize=sz2);    
    gn += 2;
  tight_layout()

######################
# test hot/cold-spot spacing (IP3R changes) - HSCSSPACEBatch - old batch - IGNORE!!!

#
def HSCSSPACEBatchLims ():
  return numpy.linspace(2,100,50)

# params for hot/cold-spot variations in IP3R and ER levels
def HSCSSPACEBatchParams (lspc=HSCSSPACEBatchLims()):
  lsec,lopt,lval = [],[],[]
  for spc in lspc: # variations in IP3R hotspot spacing
    lsec.append(["run","run","run","set","set"])
    lopt.append(["simstr","runit","saveout","ip3_origin","boost_every"])
    simstr = "HSCSSPACEBatch_HS1.7IP3R_BoostEvery_"+str(spc)
    lval.append([simstr,"1","1",str(1.7*120400.0), str(spc)])
  for spc in lspc: # variations in IP3R coldspot spacing
    lsec.append(["run","run","run","set","set"])
    lopt.append(["simstr","runit","saveout","ip3_origin","boost_every"])
    simstr = "HSCSSPACEBatch_CS0.3IP3R_BoostEvery_"+str(spc)
    lval.append([simstr,"1","1",str(0.3*120400.0), str(spc)])
  return lsec, lopt, lval

# run this batch
def HSCSSPACEBatchRun (skip=[],blog="HSCSSPACEBatch/HSCSSPACEBatch.log",bdir="HSCSSPACEBatch",qsz=12):
  batchRun(HSCSSPACEBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read data from this batch
def HSCSSPACEBatchRead ():
  return shortRead(HSCSSPACEBatchParams)

# draw data from this batch
def HSCSSPACEBatchDraw (dat=None,lnq=None,ldw=None,lsp=None,ldur=None,lamp=None,lon=None,xl=(3,102)):
  if dat is None: dat = HSCSSPACEBatchRead()
  if lnq is None: lnq,ldw,lsp,ldur,lamp,lon=getwavenqs(HSCSSPACEBatchParams,dat)
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f", "g", "h"]; lidx = [4, 20, 4+50, 20+50];
  lspc = HSCSSPACEBatchLims()
  print "a:",lidx[0],lspc[lidx[0]],"b:",lidx[1],lspc[lidx[1]]
  G = gridspec.GridSpec(2, 8);
  for idx in lidx:
    ax=subplot(G[0,gn*2:gn*2+2]); xlabel('Time(s)');
    imshow(dat[idx]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((300,700)); xlim((0,6));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
  ly = [ (4,16), (30,80), (0.9,1.2), (0.0018,0.0022) ];
  slc = [ [0,50] , [50, 101] ]; sty1 = ['r2', 'b3']; sty2 = ['rd', 'bd'];# hotspots in red, coldspots in blue
  ltitles = ["Onset (ms)", r"Speed ($\mu$m/s)", "Duration (s)", "Amplitude (mM)"];
  lletters = ["e", "f", "g", "h"]; lvals = [ lon, lsp, ldur, lamp ]; gn = 0; sz1,sz2=10,10;
  scaley = [1.0,1.0,1.0/1e3,1.0];
  for pdx in xrange(4):
    ax=subplot(G[1,gn:gn+2]); xlim(xl); xlabel("Branch Spacing (um)"); title(ltitles[pdx]); ylim(ly[pdx]);
    text(tx,ty,lletters[pdx],fontsize=14,fontweight='bold',transform=ax.transAxes); grid(True);
    for i in xrange(2):
      plot(lspc,numpy.array(lvals[pdx][slc[i][0]:slc[i][1]])*scaley[pdx],sty1[i],markersize=sz1); 
      gtt,ltt=[lspc[lidx[0]],lspc[lidx[1]]],[lvals[pdx][lidx[2*i]],lvals[pdx][lidx[2*i+1]]];
      plot(gtt,numpy.array(ltt)*scaley[pdx],sty2[i],markersize=sz2);    
    gn += 2;
  tight_layout()

######################
# test hot/cold-spot spacing (IP3R changes) - HSSPACEBatch - this is the newer batch for IP3R hotspot spacing!!!

#
def HSIP3RSPACEBatchLims ():
  return numpy.linspace(11,100,90)

# params for hot/cold-spot variations in IP3R and ER levels
def HSIP3RSPACEBatchParams (lspc=HSIP3RSPACEBatchLims()):
  lsec,lopt,lval = [],[],[]
  alls,allo,allv = [],[],[] # params shared across sims in this batch
  NewParam(alls,allo,allv,"set","ip3_stim","1.25")
  NewParam(alls,allo,allv,"set","ip3_stimT","2e3")
  NewParam(alls,allo,allv,"set","gleak","6.02")
  NewParam(alls,allo,allv,"set","gserca","2.2010625")
  NewParam(alls,allo,allv,"set","boost_halfw","5.0")
  NewParam(alls,allo,allv,"set","ip3_notorigin",str(0.5*120400.0))
  NewParam(alls,allo,allv,"set","ip3_origin",str(2.5*120400.0))
  NewParam(alls,allo,allv,"run","runit","1.0")
  NewParam(alls,allo,allv,"run","saveout","1.0")
  for spc in lspc: # variations in IP3R hotspot spacing
    tmpsec,tmpopt,tmpval=[],[],[]
    simstr = "HS_IP3R_SPACE_Batch_HS_2.5_IP3R_boost_every_"+str(spc)
    NewParam(tmpsec,tmpopt,tmpval,"run","simstr",simstr)
    for i in xrange(len(alls)): tmpsec.append(alls[i]); tmpopt.append(allo[i]); tmpval.append(allv[i]);
    NewParam(tmpsec,tmpopt,tmpval,"set","boost_every",str(spc))
    lsec.append(tmpsec); lopt.append(tmpopt); lval.append(tmpval);    
  return lsec, lopt, lval

# run this batch
def HSIP3RSPACEBatchRun (skip=[],blog="HSIP3RSPACEBatch/HSIP3RSPACEBatch.log",bdir="HSIP3RSPACEBatch",qsz=12):
  batchRun(HSIP3RSPACEBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read data from this batch
def HSIP3RSPACEBatchRead ():
  return shortRead(HSIP3RSPACEBatchParams)

# get the NQS with batch data / wave properties
def HSIP3RSPACEBatchNQS ():
  return NQS("data/13may3_ISIP3RSPACEBATCH_wp2.nqs")

# print spacing param at the specified row
def HSIP3RSPACEBatchParamsAtRow (nq,row):
  nq.tog("DB")
  SPC = nq.getcol("boost_every").x[row]
  print "row ", row, " SPACE: ", SPC

# draw the waves for first part of the figure
def HSIP3RSPACEBatchDrawWaves (nq,zoom=False):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"];
  # lidx = [14, 39, 89] 
  lidx = [9, 39, 89] 
  SPC = HSIP3RSPACEBatchLims()
  nrows = 1
  if zoom: nrows = 2
  for idx in lidx:
    ax=subplot(nrows,len(lidx),gn+1);
    dat = loadfromnq(nq,idx);
    HSIP3RSPACEBatchParamsAtRow(nq,idx) # print info
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((0,1000)); xlim((2,13));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    if not zoom: xlabel('Time (s)');
    if zoom:
      ax=subplot(nrows,len(lidx),gn+1+len(lidx)); 
      if gn == 0: ylabel(r'Position ($\mu$m)');
      xlabel('Time (s)');
      imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
      ylim((500,600)); xlim((2,5));
    gn += 1;  

# draw data from this batch - need to redo once batch data ready!!
def HSIP3RSPACEBatchDrawWaveProps (nq,xl=(3,102)):
  gn,i = 0,0; tx, ty = -0.15, 1.06
  lcols = ["speed", "dist", "time"]
  ltitles = [r'speed ($\mu$m/s)',r'distance ($\mu$m)','time (s)'];
  lspc = HSIP3RSPACEBatchLims()
  ly = [ (4,16), (30,80), (0.9,1.2), (0.0018,0.0022) ];
  sty1 = ['ko']; sty2 = ['rd', 'bd'];
  ltitles = [r"Speed ($\mu$m/s)", r"Distance ($\mu$m)", "Duration (s)", ];
  lletters = ["a", "b", "c", "d"]; gn = 0; sz1,sz2=10,10;
  scaley = [1.0,1.0,1.0/1e3]; lcol = ["speed","dist","time"]
  nq.tog("DB")
  for pdx in xrange(3):
    ax=subplot(1,3,pdx+1);
    xlim(xl); xlabel("Spacing (um)"); title(ltitles[pdx]); # ylim(ly[pdx]);
    text(tx,ty,lletters[pdx],fontsize=14,fontweight='bold',transform=ax.transAxes); grid(True);
    plot(lspc,scaley[pdx]*numpy.array(nq.getcol(lcol[pdx]).to_python()),'ko',markersize=sz1);
    plot(lspc,scaley[pdx]*numpy.array(nq.getcol(lcol[pdx]).to_python()),'k',markersize=sz1);
  # tight_layout()

######################
# gserca variations

#
def gSERCABatchLims ():
  return numpy.linspace(0.5,1.1,45)

# params for slide 24 - part that varies gserca
def gSERCABatchParams (lfctr=gSERCABatchLims()):
  lsec,lopt,lval = [],[],[]
  for fctr in lfctr:
    lsec.append(["run","run","run","set","set","set","set"])
    lopt.append(["simstr","runit","saveout","gserca","gleak","ip3_stim","ip3_stimT"])
    simstr = "gSERCABatch_"+str(fctr*1.9565)
    lval.append([simstr,"1","1",str(fctr*1.9565),"18.06","1.25","2000"])
  return lsec, lopt, lval

# run sims for gserca part of slide 24 and save output 
def gSERCABatchRun (skip=[],blog="gSERCABatch/gSERCABatch.log",bdir="gSERCABatch",qsz=12):
  batchRun(gSERCABatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read results from slide 24 sims
def gSERCABatchRead ():
  return shortRead(gSERCABatchParams)

# draw figure showing sensitivity to IP3R density
def gSERCABatchDraw (dat=None,lnq=None,ldw=None,lsp=None,ldur=None,lamp=None,lon=None,xl=(0.6,1.1),stimT=2000, lidx=[12,39]):
  if dat is None: dat = gSERCABatchRead()
  if lnq is None: lnq,ldw,lsp,ldur,lamp,lon=getwavenqs(gSERCABatchParams,dat)
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"]#; lidx = [12, 39];
  gfctr = gSERCABatchLims()
  #print "a:",lidx[0],gfctr[lidx[0]],"b:",lidx[1],gfctr[lidx[1]]
  # print the values of red dots
  naxbins = 6
  G = gridspec.GridSpec(2, 4)
  for idx in lidx:
    ax=subplot(G[0,gn*2:gn*2+2]);
    xlabel('Time(s)');
    imshow(dat[idx]["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,30,0,1000),origin='lower')
    ylim((150,850)); xlim((2,7));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
  plon = subxgtzero(lon,stimT)
  ax=subplot(G[1,0]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins); # ylim((5,20)); 
  text(tx,ty,"c",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr[12:],plon[12:],'ko',markersize=10); xlabel("SERCA density"); title("Onset (ms)"); ylim((0,350));
  gtt,ltt = [gfctr[lidx[0]], gfctr[lidx[1]] ], [plon[lidx[0]], plon[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=12);
  print '\nred dots values for onset:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nSERCAdensity={1}, onset(ms)={2}'.format(val[0],val[1],ltt[val[0]])
  ax=subplot(G[1,1]); xlim(xl);grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"d",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr[12:],lsp[12:],'ko',markersize=10); xlabel("SERCA density"); title(r"Speed ($\mu$m/s)"); ylim((0,90));
  gtt,ltt = [gfctr[lidx[0]], gfctr[lidx[1]] ], [lsp[lidx[0]], lsp[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=12);
  print '\nred dots values for speed:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nSERCAdensity={1}, speed(um/s)={2}'.format(val[0],val[1],ltt[val[0]])
  ax=subplot(G[1,2]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"e",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr[12:],numpy.array(ldur)[12:]/1e3,'ko',markersize=10); xlabel("SERCA density"); title("Duration (s)");
  gtt,ltt = [gfctr[lidx[0]], gfctr[lidx[1]] ], [ldur[lidx[0]]/1e3, ldur[lidx[1]]/1e3 ]; plot(gtt,ltt,'ro',markersize=12);
  print '\nred dots values for duration:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nSERCAdensity={1}, duration(s)={2}'.format(val[0],val[1],ltt[val[0]])
  ylim((0.0,1.45));
  ax=subplot(G[1,3]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"f",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr[12:],lamp[12:],'ko',markersize=10); xlabel("SERCA density"); title("Amplitude (mM)"); ylim((0,0.0018));
  gtt,ltt = [gfctr[lidx[0]], gfctr[lidx[1]] ], [lamp[lidx[0]], lamp[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=12);
  print '\nred dots values for amplitude:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nSERCAdensity={1}, amplitude(mM)={2}'.format(val[0],val[1],ltt[val[0]])
  tight_layout()

#####################################################################################

#####################################################################################
#                    long batch varying BOTH gserca and gip3 independently
# gip3rsercabatch

# params for sims that vary gip3 and serca
def gIP3RSERCABatchParams (lfctr1=gIP3RBatchLims(),lfctr2=gSERCABatchLims()):
  lsec,lopt,lval = [],[],[]
  for fctr1 in lfctr1[0:86]:
    for fctr2 in lfctr2[12:]:
      lsec.append(["run","run","run","run","set","set","set","set","set","set"])
      lopt.append(["simstr","runit","saveout","tstop","ip3_origin","ip3_notorigin","ip3_stim","ip3_stimT","gleak","gserca"])
      simstr = "gIP3RSERCABatch_"+str(fctr1*120400.0)+"_"+str(1.9565*fctr2)
      lval.append([simstr,"1","1","15e3",str(fctr1*120400.0),str(fctr1*120400.0),"1.25","2e3","18.06",str(1.9565*fctr2)])
  return lsec, lopt, lval

# run sims for gip3rserca and save output 
def gIP3RSERCABatchRun (skip=[],blog="gIP3RSERCABatch/gIP3RSERCABatch.log",bdir="gIP3RSERCABatch",qsz=18):
  batchRun(gIP3RSERCABatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read results from sims
def gIP3RSERCABatchRead ():
  return shortRead(gIP3RSERCABatchParams)

# makes the NQS with data from the batch - takes many hours!
def gIP3RSERCABatchMakeNQP ():
  nqp = params2nq(gIP3RSERCABatchParams)
  caCYT_init = 0.0001
  addwavenqcol(nqp,2*caCYT_init) # started at 17:17
  addwavepropcols(nqp)
  return nqp

#
def gIP3RSERCABatchSetupDat ():
  nqp = NQS("data/14sep12_gIP3RSERCABatch_wp2.nqs")
  lfctr1=gIP3RBatchLims()
  lfctr2=gSERCABatchLims()
  mspeed,mdur,mtime,mamp,monset,mdist=[numpy.zeros((86,33)) for i in xrange(6)]
  lmat = [mspeed,mdur,mtime,mamp,monset,mdist]; lcol = ['speed','dur','time','amp','onset','dist']
  for i,fctr1 in enumerate(lfctr1[0:86]):
    for j,fctr2 in enumerate(lfctr2[12:]):
      simstr = "gIP3RSERCABatch_"+str(fctr1*120400.0)+"_"+str(1.9565*fctr2)
      if nqp.select(-1,'simstr',h.SEQ,simstr) != 1: continue
      idx = int(nqp.ind.x[0])
      for mat,col in zip(lmat,lcol): mat[i,j] = nqp.getcol(col).x[idx]
  return lfctr1,lfctr2,lmat,nqp

def gIP3RSERCABatchWaves (nqp=None,lmat=None,lfctr1=gIP3RBatchLims(),lfctr2=gSERCABatchLims()):
  naxbins=6; tx, ty = -0.15, 1.06; tlx = ["a", "b", "c", "d", "e","f","g"]; txfsz=16
  if nqp is None or lmat is None: lfctr1,lfcr2,lmat,nqp = gIP3RSERCABatchSetupData()
  mspeed,mdur,mtime,mamp,monset,mdist=lmat
  # (A) top-left with 0 velocity
  # (B) top-left with > 0 velocity
  # (C) middle region
  coords = [(lfctr2[15],lfctr1[79]),(lfctr2[22],lfctr1[79]),(lfctr2[43],lfctr1[79])]
  lgs,lip3r=[],[]
  for c in coords:
    lgs.append(c[0]*1.9565)
    lip3r.append(c[1]*120400.0)
  ldat = []
  gn = 0
  for gs,ip3r in zip(lgs,lip3r):
    if nqp.select(-1,'gserca',gs,'ip3_origin',ip3r) != 1:
      print 'not found!!'
      continue
    idx = int(nqp.ind.x[0]) 
    dat = loadfromnq(nqp,idx) 
    ax=figure()
    ax.locator_params(nbins=naxbins-2);
    imshow(dat["cytca"],vmin=0,vmax=0.002,aspect='auto',extent=(0,15,0,1000),origin='lower')
    ylim((200,800)); xlim((2,8));
    #if gn == 0 or gn == 3: ylabel(r'Position ($\mu$m)');
    if gn == 0: title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=txfsz,fontweight='bold',transform=ax.transAxes)
    #if gn > 2: xlabel('Time (s)');
    ylabel(r'Position ($\mu$m)');
    gn+=1
  xlabel('Time (s)');

#
def gIP3RSERCABatchDraw (nqp=None,lmat=None,lfctr1=gIP3RBatchLims(),lfctr2=gSERCABatchLims()):
  naxbins=6; tx, ty = -0.15, 1.06; tlx = ["a", "b", "c", "d", "e","f","g"]; txfsz=16
  if nqp is None or lmat is None: lfctr1,lfcr2,lmat,nqp = gIP3RSERCABatchSetupData()
  mspeed,mdur,mtime,mamp,monset,mdist=lmat
  # (A) top-left with 0 velocity
  # (B) top-left with > 0 velocity
  # (C) middle region  
  def drawletters ():
    coords = [(lfctr2[14],lfctr1[79]),(lfctr2[17],lfctr1[69]),(lfctr2[31],lfctr1[37]),(lfctr2[40],lfctr1[8])]
    x0,x1=lfctr2[12],lfctr2[-1];
    y1,y0=lfctr1[86],lfctr1[0];
    xwid = x1-x0; ywid = y1-y0; 
    lr = ["1","2","3","4"] # text
    clrs = ["white","black","black","white"]
    for i,c in enumerate(coords):
      print i, c
      text( (c[0]-x0)/xwid, (c[1]-y0)/ywid,lr[i],color=clrs[i],\
          horizontalalignment='center',verticalalignment='center',fontsize=18,fontweight='bold',transform=ax.transAxes);
  gn = 0
  ax=subplot(2,2,1); ax.locator_params(nbins=naxbins);
  title("Onset (ms)"); ylabel(r'$IP_3R$ density'); drawletters()
  imshow(monset-2e3,vmin=0,vmax=400,origin='lower',extent=(lfctr2[12],lfctr2[-1],lfctr1[0],lfctr1[86]),aspect='auto',interpolation='None'); colorbar(); text(tx,ty,tlx[gn],fontsize=txfsz,fontweight='bold',transform=ax.transAxes); gn+=1
  ax=subplot(2,2,2); ax.locator_params(nbins=naxbins);
  imshow(mspeed,origin='lower',extent=(lfctr2[12],lfctr2[-1],lfctr1[0],lfctr1[86]),aspect='auto',interpolation='None'); colorbar()
  title(r'Speed ($\mu$m/s)'); drawletters()
  text(tx,ty,tlx[gn],fontsize=txfsz,fontweight='bold',transform=ax.transAxes); 
  gn+=1
  ax=subplot(2,2,3); ax.locator_params(nbins=naxbins);
  imshow(mdur/1e3,vmin=0,vmax=4,origin='lower',extent=(lfctr2[12],lfctr2[-1],lfctr1[0],lfctr1[86]),aspect='auto',interpolation='None'); colorbar()
  xlabel('SERCA density'); ylabel(r'$IP_3R$ density'); title("Duration (s)"); drawletters()
  text(tx,ty,tlx[gn],fontsize=txfsz,fontweight='bold',transform=ax.transAxes); gn+=1
  ax=subplot(2,2,4); ax.locator_params(nbins=naxbins);
  imshow(mamp,vmin=.0014,vmax=.0017,origin='lower',extent=(lfctr2[12],lfctr2[-1],lfctr1[0],lfctr1[86]),aspect='auto',interpolation='None'); colorbar()
  xlabel('SERCA density'); title("Amplitude (mM)"); drawletters()
  text(tx,ty,tlx[gn],fontsize=txfsz,fontweight='bold',transform=ax.transAxes); gn+=1

#####################################################################################
# AMPA stim (ER priming) variations

#
def AMPAStimBatchLims (start=0, stop=150, size=31):
  a = list(numpy.linspace(start,stop,size))
  #a.append(1.0) # commented this out to run batches to plot fig 11
  a.sort()
  return numpy.array(a)

# params for part that varies AMPA stim
def AMPAStimBatchParams (lnsn=AMPAStimBatchLims(0,150,2)):
  lsec,lopt,lval = [],[],[]
  for nstimNumber in lnsn:
    lsec.append(["run","run","run","run","set","set","set","set","set","set","set","set","set","set","run","run","run","run"])
    lopt.append(["simstr","runit","cvodeactive","saveout","nstimNumber","ip3_stim","ip3_stimT","ip3_init","gCaChannels","nstimStart","nstimInterval","nstimNumber","nconnWeight","nconnActive","loadstate","statestr","loadRXDStateOnly","tstop"])
    simstr = "AMPAStimBatch_"+str(int(nstimNumber))
    lval.append([simstr,"1","1","1",str(int(nstimNumber)),"2.5","7e3","0.0","1.e-6","3e3","25",str(nstimNumber),"0.5","1","1","14aug12_3000s_save_C","1","12000"])
  return lsec, lopt, lval

# run sims 
def AMPAStimBatchRun (skip=[],blog="AMPAStimBatch/AMPAStimBatch.log",bdir="AMPAStimBatch",qsz=20):
  batchRun(AMPAStimBatchParams,blog,skip=skip,qsz=qsz,bdir=bdir)

# read results 
def AMPAStimBatchRead ():
  return shortRead(AMPAStimBatchParams)

# draw figure showing sensitivity to number of AMPA stims
def AMPAStimBatchDrawOLD (dat=None,lnq=None,ldw=None,lsp=None,ldur=None,lamp=None,lon=None,xl=(0,150),stimT=7000, lidx=[1,31]):
  if dat is None: dat = AMPAStimBatchRead()
  if lnq is None: lnq,ldw,lsp,ldur,lamp,lon=getwavenqs(AMPAStimBatchParams,dat)
  gn,i = 0,0; tx, ty = -0.15, 1.06
  tlx = ["a", "b", "c", "d", "e", "f"]#; lidx = [12, 39];
  gfctr = AMPAStimBatchLims()
  #print "a:",lidx[0],gfctr[lidx[0]],"b:",lidx[1],gfctr[lidx[1]]
  # print the values of red dots
  naxbins = 6
  G = gridspec.GridSpec(2, 4)
  for idx in lidx:
    ax=subplot(G[0,gn*2:gn*2+2]);
    xlabel('Time(s)');
    imshow(dat[idx]["cytca"],vmin=0,vmax=0.0075,aspect='auto',extent=(0,12,0,1000),origin='lower')
    ylim((375,625)); xlim((6.5,10));
    if gn == 0: ylabel(r'Position ($\mu$m)');
    title(r'$Ca^{2+}_{cyt}$ (mM)'); 
    text(tx,ty,tlx[gn],fontsize=14,fontweight='bold',transform=ax.transAxes)
    gn += 1;
  plon = subxgtzero(lon,stimT)
  ax=subplot(G[1,0]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins); # ylim((5,20)); 
  text(tx,ty,"c",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,plon,'ko',markersize=10); xlabel("AMPA stims"); title("Onset (ms)"); #ylim((0,350));
  gtt,ltt = [gfctr[lidx[0]], gfctr[lidx[1]] ], [plon[lidx[0]], plon[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=12);
  print '\nred dots values for onset:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nAMPAStims={1}, onset(ms)={2}'.format(val[0],val[1],ltt[val[0]])
  ax=subplot(G[1,1]); xlim(xl);grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"d",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,lsp,'ko',markersize=10); xlabel("AMPA stims"); title(r"Speed ($\mu$m/s)"); #ylim((0,90));
  gtt,ltt = [gfctr[lidx[0]], gfctr[lidx[1]] ], [lsp[lidx[0]], lsp[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=12);
  print '\nred dots values for speed:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nAMPAStims={1}, speed(um/s)={2}'.format(val[0],val[1],ltt[val[0]])
  ax=subplot(G[1,2]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"e",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,numpy.array(ldur)/1e3,'ko',markersize=10); xlabel("AMPA stims"); title("Duration (s)");
  gtt,ltt = [gfctr[lidx[0]], gfctr[lidx[1]] ], [ldur[lidx[0]]/1e3, ldur[lidx[1]]/1e3 ]; plot(gtt,ltt,'ro',markersize=12);
  print '\nred dots values for duration:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nAMPAStims={1}, duration(s)={2}'.format(val[0],val[1],ltt[val[0]])
  #ylim((0.0,1.45));
  ax=subplot(G[1,3]); xlim(xl); grid(True); ax.locator_params(nbins=naxbins);
  text(tx,ty,"f",fontsize=14,fontweight='bold',transform=ax.transAxes);
  plot(gfctr,lamp,'ko',markersize=10); xlabel("AMPA stims"); title("Amplitude (mM)"); #ylim((0,0.0018));
  gtt,ltt = [gfctr[lidx[0]], gfctr[lidx[1]] ], [lamp[lidx[0]], lamp[lidx[1]] ]; plot(gtt,ltt,'ro',markersize=12);
  print '\nred dots values for amplitude:'
  for val in enumerate(gtt):
    print 'red dot index({0}):\nAMPAStims={1}, amplitude(mM)={2}'.format(val[0],val[1],ltt[val[0]])
  tight_layout()

#
def mynorm (x):
  y = x - amin(x)
  y = y / amax(y)
  return y

# draw figure showing sensitivity to number of AMPA stims
def AMPAStimBatchDraw (dat=None,lnq=None,ldw=None,lsp=None,ldur=None,lamp=None,lon=None,xl=(2.5,10),yl=(400,600),lidx=[0,31]):
  if dat is None: dat = AMPAStimBatchRead()
  #if lnq is None: lnq,ldw,lsp,ldur,lamp,lon=getwavenqs(AMPAStimBatchParams,dat)
  labsz=18; titlesz=18;
  gn,i = 1,0; tx, ty = -0.15, 1.06
  nrow,ncol=3,3; 
  gfctr = AMPAStimBatchLims(); naxbins = 6
  ld = dat[lidx[0]]
  if ldur is not None:
    print lidx[0],'input:',gfctr[lidx[0]],'speed:',lsp[lidx[0]], 'dur:',ldur[lidx[0]]/1e3, 'amp:',lamp[lidx[0]], 'onset:',lon[lidx[0]]
  ax=subplot(nrow,ncol,1);
  text(tx,ty,"a",fontsize=labsz,fontweight='bold',transform=ax.transAxes);
  imshow(ld['volt'],vmin=-76.5,vmax=-8,aspect='auto',origin='lower',extent=(0,12,0,1000)); colorbar()
  title(r'Voltage (mV)',fontsize=titlesz); 
  ylabel(r'Position ($\mu$m)',fontsize=labsz);
  ax=subplot(nrow,ncol,4);
  imshow(ld['erca'],vmin=0.006,vmax=0.06,aspect='auto',origin='lower',extent=(0,12,0,1000)); colorbar()
  ylabel(r'Position ($\mu$m)',fontsize=labsz);
  title(r'$Ca^{2+}_{ER}$ (mM)',fontsize=titlesz); 
  ax=subplot(nrow,ncol,7); 
  imshow(ld['cytca'],vmin=0.00005,vmax=0.0075,aspect='auto',origin='lower',extent=(0,12,0,1000)); colorbar()
  title(r'$Ca^{2+}_{cyt}$ (mM)',fontsize=titlesz); 
  xlabel('Time(s)',fontsize=labsz);
  ylabel(r'Position ($\mu$m)',fontsize=labsz);
  ld = dat[lidx[1]]
  if ldur is not None:
    print lidx[1],'input:',gfctr[lidx[1]],'speed:',lsp[lidx[1]], 'dur:',ldur[lidx[1]]/1e3, 'amp:',lamp[lidx[1]], 'onset:',lon[lidx[1]]
  ax=subplot(nrow,ncol,2);
  text(tx,ty,"b",fontsize=labsz,fontweight='bold',transform=ax.transAxes);
  imshow(ld['volt'],vmin=-76.5,vmax=-8,aspect='auto',origin='lower',extent=(0,12,0,1000)); colorbar()
  title(r'Voltage (mV)',fontsize=titlesz); 
  ax=subplot(nrow,ncol,5);
  imshow(ld['erca'],vmin=0.006,vmax=0.06,aspect='auto',origin='lower',extent=(0,12,0,1000)); colorbar()
  title(r'$Ca^{2+}_{ER}$ (mM)',fontsize=titlesz); 
  ax=subplot(nrow,ncol,8);
  title(r'$Ca^{2+}_{cyt}$ (mM)',fontsize=titlesz); 
  imshow(ld['cytca'],vmin=0.00005,vmax=0.0075,aspect='auto',origin='lower',extent=(0,12,0,1000)); colorbar()
  xlabel('Time(s)',fontsize=labsz);
  for i in [1,2,4,5,7,8]:
    ax=subplot(nrow,ncol,i); ax.locator_params(nbins=naxbins);
    xlim(xl); ylim(yl); 
  ta = np.arange(0,12e3,5) / 1e3;
  ax=subplot(nrow,ncol,3); ax.locator_params(nbins=naxbins);
  text(tx,ty,"c",fontsize=labsz,fontweight='bold',transform=ax.transAxes);
  plot(ta,dat[lidx[0]]['volt'][499,:],'k')
  plot(ta,dat[lidx[1]]['volt'][499,:],'r',linewidth=1)
  title(r'Voltage (mV)',fontsize=titlesz); 
  xlim(xl)
  ax=subplot(nrow,ncol,6); ax.locator_params(nbins=naxbins);
  plot(ta,dat[lidx[0]]['erca'][499,:],'k')
  plot(ta,dat[lidx[1]]['erca'][499,:],'r',linewidth=1)
  title(r'$Ca^{2+}_{ER}$ (mM)',fontsize=titlesz); 
  xlim(xl)
  ax=subplot(nrow,ncol,9); ax.locator_params(nbins=naxbins);
  plot(ta,dat[lidx[0]]['cytca'][499,:],'k')
  plot(ta,dat[lidx[1]]['cytca'][499,:],'r',linewidth=1)
  title(r'$Ca^{2+}_{cyt}$ (mM)',fontsize=titlesz); 
  xlim(xl); ylim((0.00005,0.00015));
  xlabel('Time(s)',fontsize=labsz);

def AMPAStimBatchDraw_fig11 (dat=None,xl=(2.5,10),yl=(400,600)):#,lidx=[0,31]):
  if dat is None: dat = AMPAStimBatchRead()
  #if lnq is None: lnq,ldw,lsp,ldur,lamp,lon=getwavenqs(AMPAStimBatchParams,dat)
  gn,i = 1,0; tx, ty = -0.15, 1.06
  nrow,ncol=3,3;  ex = (0,12,0,1000); naxbins = 6
  spid = 0 #subplot id
  for ld in dat:
#      ld = dat[lidx[0]];  
#      ax=subplot(nrow,ncol,1);
      ax=subplot(nrow,ncol,spid+1);
      text(tx,ty,"a",fontsize=14,fontweight='bold',transform=ax.transAxes);
      imshow(ld['volt'],vmin=-76.5,vmax=-8,aspect='auto',origin='lower',extent=ex); colorbar()
      title(r'Voltage (mV)'); 
      ylabel(r'Position ($\mu$m)');
#      ax=subplot(nrow,ncol,4);
      ax=subplot(nrow,ncol,spid+4);      
      imshow(ld['erca'],vmin=0.006,vmax=0.06,aspect='auto',origin='lower',extent=ex); colorbar()
      ylabel(r'Position ($\mu$m)');
      title(r'$Ca^{2+}_{ER}$ (mM)'); 
#      ax=subplot(nrow,ncol,7); 
      ax=subplot(nrow,ncol,spid+7); 
      imshow(ld['cytca'],vmin=0.00005,vmax=0.0075,aspect='auto',origin='lower',extent=ex); colorbar()
      title(r'$Ca^{2+}_{cyt}$ (mM)'); 
      xlabel('Time(s)');
      ylabel(r'Position ($\mu$m)');
      spid +=1
      # ld = dat[lidx[1]]
      # ax=subplot(nrow,ncol,2);
      # text(tx,ty,"b",fontsize=14,fontweight='bold',transform=ax.transAxes);
      # imshow(ld['volt'],vmin=-76.5,vmax=-8,aspect='auto',origin='lower',extent=ex); colorbar()
      # title(r'Voltage (mV)'); 
      # ax=subplot(nrow,ncol,5);
      # imshow(ld['erca'],vmin=0.006,vmax=0.06,aspect='auto',origin='lower',extent=ex); colorbar()
      # title(r'$Ca^{2+}_{ER}$ (mM)'); 
      # ax=subplot(nrow,ncol,8);
      # title(r'$Ca^{2+}_{cyt}$ (mM)'); 
      # imshow(ld['cytca'],vmin=0.00005,vmax=0.0075,aspect='auto',origin='lower',extent=ex); colorbar()
      # xlabel('Time(s)');
  for i in [1,2,4,5,7,8]:
    ax=subplot(nrow,ncol,i); ax.locator_params(nbins=naxbins);
    xlim(xl); ylim(yl); 
  ta = np.arange(0,12e3,5) / 1e3;
  ax=subplot(nrow,ncol,3); ax.locator_params(nbins=naxbins);
  text(tx,ty,"c",fontsize=14,fontweight='bold',transform=ax.transAxes);
  for pair in zip(dat, ['k', 'r']):
    plot(ta,pair[0]['volt'][499,:],pair[1])
  # plot(ta,dat[lidx[0]]['volt'][499,:],'k')
  # plot(ta,dat[lidx[1]]['volt'][499,:],'r',linewidth=1)
  title(r'Voltage (mV)'); 
  xlim(xl)
  ax=subplot(nrow,ncol,6); ax.locator_params(nbins=naxbins);
  for pair in zip(dat, ['k', 'r']):
    plot(ta,pair[0]['erca'][499,:],pair[1])
  # plot(ta,dat[lidx[0]]['erca'][499,:],'k')
  # plot(ta,dat[lidx[1]]['erca'][499,:],'r',linewidth=1)
  title(r'$Ca^{2+}_{ER}$ (mM)'); 
  xlim(xl)
  ax=subplot(nrow,ncol,9); ax.locator_params(nbins=naxbins);
  for pair in zip(dat, ['k', 'r']):
    plot(ta,pair[0]['cytca'][499,:],pair[1])
  # plot(ta,dat[lidx[0]]['cytca'][499,:],'k')
  # plot(ta,dat[lidx[1]]['cytca'][499,:],'r',linewidth=1)
  title(r'$Ca^{2+}_{cyt}$ (mM)'); 
  xlim(xl); ylim((0.00005,0.00015));
  xlabel('Time(s)');




#####################################################################################

# get wave features using dat and lnq
def getwavepropl (dat,lnq,nooverlap=True,uponly=True):
  ldw,lsp,ldur,lamp,lon = [],[],[],[],[] # list of dictionaries, speed, duration, amplitude, onsets
  i = 0
  for nq in lnq:
    nq.verbose=0
    nq.select("widx",1,"y",">=",500)
    if nooverlap: nq.select("&&","overlap",0)
    if uponly: nq.select("&&","up",1)
    ldw.append(nqs2pyd(nq)); 
    lsp.append(wavespeed(dat[i],key="cytca",nqw=nq))
    ldur.append(wavedur(dat[i],key="cytca",nqw=nq))
    lamp.append(amax(dat[i]["cytca"]))
    lon.append(waveonset(dat[i],key="cytca",nqw=nq))
    nq.verbose=1
    i += 1
  return ldw,lsp,ldur,lamp,lon

# read the nqs objects made with getwavenqs previously saved
def readlnq (whichParams,dat=None):
  lsec,lopt,lval = whichParams()
  if dat is None: dat = shortRead(whichParams)
  simstridx = lopt[0].index('simstr')
  lnq = []
  for i in xrange(len(dat)):
    fn = "./data/" + lval[i][simstridx] + "_wnq.nqs"
    nq = NQS(fn)
    lnq.append(nq)
  return lnq

# save the nqs objects made with getwavenqs
def savelnq (whichParams,lnq):
  lsec,lopt,lval = whichParams()
  i=0
  simstridx = lopt[0].index('simstr')
  for nq in lnq: 
    nq.tog("DB")
    fn = "./data/" + lval[i][simstridx] + "_wnq.nqs" # assumes lval[i][0] contains simstr
    nq.sv(fn)
    i += 1

# get wavenqs associated with a batch
#  NB: novup specifies whether to exclude the overlapping runs and only look at
#      upper portion of the wave. this is important for good stats on wave features
#      and fine since waves typically symmetric.
def getwavenqs (whichParams,thresh=caCYT_init*2,dat=None,nooverlap=True,uponly=True):
  lsec,lopt,lval = whichParams()
  if dat is None: dat = shortRead(whichParams)
  from cawave import recdt,wavenq,tstop
  lnq,ldw,lsp,ldur,lamp,lon = [],[],[],[],[],[] # list of nqs, dictionaries, speed, duration, amplitude, onsets
  for i in xrange(len(dat)):
    print i
    nq = wavenq(dat[i]["cytca"],thresh=thresh)
    nq.select("widx",1,"y",">=",500)
    if nooverlap: nq.select("&&","overlap",0)
    if uponly: nq.select("&&","up",1)
    ldw.append(nqs2pyd(nq)); lnq.append(nq);
    lsp.append(wavespeed(dat[i],key="cytca",nqw=nq))
    ldur.append(wavedur(dat[i],key="cytca",nqw=nq))
    lamp.append(amax(dat[i]["cytca"]))
    lon.append(waveonset(dat[i],key="cytca",nqw=nq))
  return lnq,ldw,lsp,ldur,lamp,lon

# ldw is a list of wave nqs from getwavenqs
def plotwavenqs (ldw,cmap=cm.Greys,lylims=[(500,750),(0,5)]):
  csm=cm.ScalarMappable(cmap=cmap); csm.set_clim((0,1));
  lwhich = ["pos", "dur"]; tx, ty = -0.15, 1.06; tlx = ["a", "b", "c"]; idx=0;
  for which in lwhich:
    ax=subplot(1,len(lwhich),idx+1)
    if which == "pos":
      for i in xrange(len(ldw)): plot(ldw[i]['startt']/1e3,ldw[i]['y'],color=csm.to_rgba(float(i)/(len(ldw))))
      for i in xrange(len(ldw)): plot(ldw[i]['startt']/1e3,ldw[i]['y'],'o',color=csm.to_rgba(float(i)/(len(ldw))))
      ylabel('Position (um)'); 
      xlim((0,25)); ylim(lylims[idx]); xlabel('Time(s)'); 
    elif which == "speed":
      pks = numpy.zeros((len(ldw),1))
      for j in xrange(len(ldw)):
        if len(ldw[j]['speed']) > 0:
          pks[j][0] = amax( ldw[j]['speed'] )
      # for i in xrange(len(ldw)): plot(ldw[i]['startt']/1e3,ldw[i]['speed'],color=csm.to_rgba(float(i)/(len(ldw))))
      # for i in xrange(len(ldw)): plot(ldw[i]['startt']/1e3,ldw[i]['speed'],'o',color=csm.to_rgba(float(i)/(len(ldw))))
      xx = numpy.linspace(0,2,len(ldw))
      plot(xx,pks,'ko'); ylabel('Peak Speed (um/ms)'); xlabel('g');
      # ylabel('Speed (um/ms)'); 
    elif which == "dur":
      for i in xrange(len(ldw)): plot(ldw[i]['startt']/1e3,ldw[i]['durt']/1e3,color=csm.to_rgba(float(i)/(len(ldw))))
      for i in xrange(len(ldw)): plot(ldw[i]['startt']/1e3,ldw[i]['durt']/1e3,'o',color=csm.to_rgba(float(i)/(len(ldw))))
      ylabel('Duration (s)');
      xlim((0,25)); ylim(lylims[idx]); xlabel('Time(s)');
    text(tx,ty,tlx[idx],fontsize=14,fontweight='bold',transform=ax.transAxes)
    idx += 1

