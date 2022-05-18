from pylab import *
from matplotlib import pyplot as mp
mp.ion() # for interactive plotting
from scipy import ndimage
from scipy.interpolate import interp1d
import sys,os,numpy,subprocess
from neuron import h,rxd
from math import ceil
h.load_file("stdrun.hoc") # creates cvode object
from nqs import *
from conf import *
from time import time, clock

from time import clock # to time simulation
import datetime # to format time of run
import pickle
import ConfigParser
dataDir = './data/'


# determine config file name
def setfcfg ():
  '''determine config file name'''
  fcfg = "cawave.cfg" # default config file name
  for i in xrange(len(sys.argv)):
    if sys.argv[i].endswith(".cfg") and os.path.exists(sys.argv[i]):
      fcfg = sys.argv[i]
  print "config file is " , fcfg
  return fcfg

fcfg=setfcfg() # config file name

# read config file, set params
er_scale,ip3_notorigin,ip3_origin,gserca0,boost_every,\
ip3_stim,gleak0,ip3rtau,recdt,tstop,cvodeactive,dt,runit,\
simstr,saveout,stim_minx,stim_maxx,dodraw,loadState,ip3_stimT,\
ip3_init,caCYT_init,caAvg_init,caDiff,ip3Diff,gCaChannels,spaceum,\
ca_stim,ca_stimT,boost_halfw,synLocl,nstimStart,nstimInterval,nstimNumber,\
nconnThreshold, nconnDelay, nconnWeight, nconnActive,\
ipydebug, statestr, saveState, useInitDict,\
IP3ForceInit, loadRXDStateOnly, electrical = readconf(fcfg) 

simdur = tstop # duration of simulation
tstart = 0
h.tstop = tstop

if electrical: print 'running with electricity elements and a single AMPA synapse'
else: print 'running without electricity'

def pnone (): print ''
debug_here = pnone

if ipydebug:
  from IPython.core.debugger import Tracer;
  debug_here = Tracer() #debugging using ipdb

def printt (s=''): print s, ': t = ',h.t

# backup the config file
def backupcfg (simstr):
  if not os.path.exists('backupcfg'): os.mkdir('backupcfg')
  fout = 'backupcfg/' + simstr + '.cfg'
  if os.path.exists(fout):
    print 'removing prior cfg file' , fout
    os.system('rm ' + fout)  
  os.system('cp ' + fcfg + ' ' + fout) # fcfg created in geom.py via conf.py

backupcfg(simstr) # backup the config file

# ip3_notorigin: initial concentratoin of IP3
# ip3_origin: concentration for non-homogenous distribution of IP3R

if cvodeactive != 0:
  h.cvode.active(1) # turn on variable time-step integration
  h.cvode.atol(1e-6)
else:
  h.cvode.active(0) # turn off variable time-step integration

# make the 'cell' - just a dend
def makecell ():
  '''make the 'cell' - just a dend - parameters obtained from Safiulina et al, 2010'''
  dend = h.Section()
  dend.diam = 1
  dend.L = 1000 # 251 um - use odd number of segmnets
  dend.nseg = int(dend.L / spaceum)
  if electrical:
    dend.Ra = 150
    Rm =  25370
    for myseg in dend: myseg.v = -72.205747407792941#-64 
    for myseg in dend: myseg.cm = 1.41
    dend.insert('pas')
    for myseg in dend: myseg.pas.g = 1.0/Rm
    for myseg in dend: myseg.pas.e = -64
    gc = gCaChannels
    dend.insert('cal') # insert L-type Ca channel
    for myseg in dend: myseg.cal.gcalbar = gc #1e-30 # 1.e-30 # gc #1.e-3 #1.e-5
    dend.insert('cagk')
    gKc = 9.e-3
    for myseg in dend: myseg.cagk.gbar=gKc
    dend.insert('can')
    for myseg in dend: myseg.can.gcanbar=gc  
    dend.insert('cat')
    for myseg in dend: myseg.cat.gcatbar = gc  
    dend.insert('kap')
    KMULTP = 0.02
    for myseg in dend: myseg.kap.gkabar = KMULTP
    dend.insert('kdr')
    for myseg in dend: myseg.kdr.gkdrbar = 0.01 
    dend.insert('na3')
    AXONM = 5
    for myseg in dend: myseg.na3.gbar = 0.022*AXONM
  return dend

h.celsius = 37
dend = makecell()

# some parameters from develop.hoc from Safiulina et al, 2010
def insertSyn(synLocl, nstimStart, nstimInterval, nstimNumber, nconnThreshold, nconnDelay, nconnWeight):
  '''will insert ExpSyn synapses at locations in synLocl, and return lists of netStim, netConn objects & of vectors recording the spike times from netStims'''
  # add list of ExpSyn synapses
  synlist = []
  #for i in arange(0, 1, 0.1):
  for i in synLocl:
    mysyn = h.Exp2Syn(i,sec=dend)
    #mysyn = h.MyExp2SynBB(i,sec=dend)
    mysyn.tau1 = 0.05
    mysyn.tau2 = 5.3
    mysyn.e = 0
    synlist.append(mysyn)  
  lrandomSeed = [1234, 2345, 3456, 4567, 5678, 6789, 7890, 8901, 9012, 9876]  
  # add list of netstim
  netStimlist = []
  for i in xrange(len(synlist)):
    myns = h.NetStim()
    myns.seed(1234)
    myns.start = nstimStart # to coincide with the IP3 stimulus
    myns.interval = nstimInterval
    myns.number = nstimNumber
    netStimlist.append(myns)  
  # add list of netcon
  netConlist = []
  netConRecveclist = [] # list of vectors which record netCon event times
  for i in xrange(len(synlist)):
    mync = h.NetCon(netStimlist[0], synlist[i], nconnThreshold, nconnDelay, nconnWeight) # one netstim conneted to all synapses 
    # zero delay (to make it instantaneous)
    mync.active(nconnActive)
    netConlist.append(mync)
    myncvec = h.Vector()
    mync.record(myncvec)# record time events from netCon
    netConRecveclist.append(myncvec)
  if len(synlist) ==0: print 'no synapses inserted'
  else: print 'inserted synapse(s) at dend location(s):', synLocl
  return synlist, netStimlist, netConlist, netConRecveclist

synlist, netStimlist, netConlist, netConRecveclist = None,None,None,None
if electrical:
  synlist, netStimlist, netConlist, netConRecveclist = insertSyn(synLocl,nstimStart,nstimInterval,nstimNumber,\
                                                                 nconnThreshold,nconnDelay,nconnWeight)

# set synapse properties based on params in config file
def setSynProps ():
  if not electrical: return
  for ns in netStimlist:
    ns.number = nstimNumber
    ns.interval = nstimInterval
    ns.start = nstimStart
  for nc in netConlist:
    nc.active(nconnActive)
    nc.weight[0] = nconnWeight
    nc.threshold = nconnThreshold
    nc.delay = nconnDelay

fc = 0.83 # fraction of cytoplasmic volume
fe = 0.17 # fraction of ER volume

# the surface_fraction=1 means that even though cyt only occupies fc of the space,
# it has all of the surface area (important for currents from NMODL)
cyt = rxd.Region(h.allsec(), nrn_region='i', geometry=rxd.FractionalVolume(fc, surface_fraction=1))
er = rxd.Region(h.allsec(), geometry=rxd.FractionalVolume(fe))

# this defines a boundary with an area independent of geometry (here there is
# one square micron of area per micron of length)
cyt_er_membrane = rxd.Region(h.allsec(), geometry=rxd.ScalableBorder(1))

# reaction-diffusion definition
ca = rxd.Species([cyt, er], d=caDiff, name='ca', charge=2, initial=caCYT_init)
ip3 = rxd.Species(cyt, d=ip3Diff, initial=ip3_init)

#gip3r - not a rate constant but multiplies k, which is a rate constant (?)
gip3r = rxd.Parameter(cyt, initial=ip3_notorigin) # parameters but set magnitude via concentration or value
gserca = rxd.Parameter(cyt, initial=gserca0) 
gleak = rxd.Parameter(cyt, initial=gleak0) 

# action of IP3 receptor
Kip3 = 0.13 # Kip3 = 0.15
Kact = 0.4
minf = ip3[cyt] * 1000. * ca[cyt] / (ip3[cyt] + Kip3) / (1000. * ca[cyt] + Kact)
#ip3r_gate_state = rxd.State(cyt_er_membrane, initial=0.78282068629348611)
ip3r_gate_state = rxd.State(cyt_er_membrane, initial=0.8)
h_gate = ip3r_gate_state[cyt_er_membrane]
k = gip3r[cyt] * (minf * h_gate) ** 3 
ip3r = rxd.MultiCompartmentReaction(ca[er]<>ca[cyt], k, k, membrane=cyt_er_membrane)

# IP3 receptor gating
ip3rg = rxd.Rate(h_gate, (1. / (1 + 1000. * ca[cyt] / (0.4)) - h_gate) / ip3rtau)

# SERCA pump: pumps ca from cyt -> ER
Kserca = 0.1 # Michaelis constant for SERCA pump
serca = rxd.MultiCompartmentReaction(ca[cyt]>ca[er],gserca[cyt]*(1e3*ca[cyt])**2/(Kserca**2+(1e3*ca[cyt])**2),membrane=cyt_er_membrane,custom_dynamics=True)

# leak channel: bidirectional ca flow btwn cyt <> ER
leak = rxd.MultiCompartmentReaction(ca[er]<>ca[cyt], gleak[cyt], gleak[cyt], membrane=cyt_er_membrane)


def getInitDict(fname=fcfg):
  '''return a dictioinary containing the state variables and their values, obtained from initation file fname'''
  config = ConfigParser.ConfigParser()
  config.read(fname)
  varnameList = ['o_cagk', 'm_cal', 'h_can', 'm_can', 'h_cat', 'm_cat', 'l_kap', 'n_kap', 'n_kdr', 'h_na3', 'm_na3', 's_na3', 'ip3', 'ip3r_gate_state', 'v', 'ca_cyt']#, 'ca_er']
  initDict = {}
  for varname in varnameList: initDict[varname] = config.getfloat('sstvars',varname)
  return initDict

## read values from config file
#initDict = getInitDict(fcfg)

def printGetEventQueueinfo():
  '''will print the size of tvec, flagvec and targetlist obtained from h.cvode.event_queue_info'''
  tvec = h.Vector()
  flagvec = h.Vector()
  targetlist = h.List()
  h.cvode.event_queue_info(3, tvec, flagvec, targetlist)
  print 'tvec:', tvec.printf()
  print 'flagvec:', flagvec.printf()
  print 'targetlist count:', targetlist.count()
  return tvec, flagvec, targetlist

def setHocStateVars(fname=fcfg):
  '''will initiatlize hoc state variables, using values from dictionary initDict'''
  initDict = getInitDict(fname)
  dend.v = initDict['v']
  dend.o_cagk = initDict['o_cagk']
  dend.m_cal = initDict['m_cal']
  dend.h_can = initDict['h_can']
  dend.m_can = initDict['m_can']
  dend.h_cat = initDict['h_cat']
  dend.m_cat = initDict['m_cat']
  dend.l_kap = initDict['l_kap']
  dend.n_kap = initDict['n_kap']
  dend.n_kdr = initDict['n_kdr']
  dend.h_na3 = initDict['h_na3']
  dend.m_na3 = initDict['m_na3']
  dend.s_na3 = initDict['s_na3']
  print 'from setHocStateVars:\n', printGetEventQueueinfo()

def setRxDStateVars(fname=fcfg):
  '''will initilize rxd state variables using values from dictionary initDict'''
  initDict = getInitDict(fname)
  caer_init = (caAvg_init - initDict['ca_cyt'] * fc) / fe
  ca[er].concentration = caer_init
  ca[cyt].concentration = initDict['ca_cyt']
  ip3[cyt].concentration = initDict['ip3']
  ip3r_gate_state[cyt_er_membrane].concentration = initDict['ip3r_gate_state']
  print 'from setRxDStateVars:\n', printGetEventQueueinfo()
  
def getStateVarsDict():
    '''will return a dictionary with hoc state variables and their values at dendloc.'''
    #myseg = dend(dendloc)
    stvardict = {}
    stvardict['cal_m'] = np.array([myseg.cal.m for myseg in dend])
    stvardict['can_m'] = np.array([myseg.can.m for myseg in dend])
    stvardict['can_h'] = np.array([myseg.can.h for myseg in dend])    
    stvardict['cat_m'] = np.array([myseg.cat.m for myseg in dend])    
    stvardict['cat_h'] = np.array([myseg.cat.h for myseg in dend])
    stvardict['kap_n'] = np.array([myseg.kap.n for myseg in dend])    
    stvardict['kap_l'] = np.array([myseg.kap.l for myseg in dend])
    stvardict['kdr_n'] = np.array([myseg.kdr.n for myseg in dend])    
    stvardict['cagk_o'] = np.array([myseg.cagk.o for myseg in dend])
    stvardict['na3_m'] = np.array([myseg.na3.m for myseg in dend])    
    stvardict['na3_h'] = np.array([myseg.na3.h for myseg in dend])    
    stvardict['na3_s'] = np.array([myseg.na3.s for myseg in dend])
    stvardict['v'] = np.array([myseg.v for myseg in dend])    
    stvardict['cai'] = np.array(ca[cyt].concentration)
    stvardict['ip3'] = np.array(ip3[cyt].concentration)
    stvardict['h_gate'] = np.array(h_gate.concentration)
    stvardict['ca_er'] = np.array(ca[er].concentration)
    return stvardict

# place the ip3 to stimulate ca wave
def place_ip3_stim (minx=stim_minx,maxx=stim_maxx,IP3Stim=ip3_stim):
  if IP3Stim <= 0.0: return # no null stim - use baseline for empty stim
  print "\ntime : " , h.t, " , placing ip3 stim from : " , minx, " to " , maxx, " . val = " , IP3Stim
  # this gives us four nodes (498.5, 499.5, 500.5, 501.5)
  for node in ip3.nodes:
    if minx < node.x * dend.L < maxx:
      node.concentration = IP3Stim
  if ip3_stimT > 0: # only re_init if called from event queue
    # rxd.re_init() # 
    h.cvode.re_init()

# check the ip3 stim for ca wave
def check_ip3_stim (minx=stim_minx,maxx=stim_maxx,IP3Stim=ip3_stim):
  print "\ntime : " , h.t, " , checking ip3 stim from : " , minx, " to " , maxx, " . val = " , IP3Stim
  # this gives us four nodes (498.5, 499.5, 500.5, 501.5)
  stimxlist = []
  for node in ip3.nodes:
    if minx < node.x * dend.L < maxx:
      print "node.x " , node.x, ", concentration is " , node.concentration
      stimxlist.append(node.x)

def get_ip3stimLoc(minx=stim_minx, maxx=stim_maxx):
  '''return locations on dendrite where ip3stim has been put'''
  stimxlist = []
  for node in ip3.nodes:
    if minx < node.x * dend.L <maxx:
      stimxlist.append(node.x)
  return stimxlist

# place ca to stimulate ca wave
def place_ca_stim (minx=stim_minx,maxx=stim_maxx,CAStim=ca_stim):
  if CAStim <= 0.0: return # no null stim - use baseline for empty stim
  print "\ntime : " , h.t, " , placing ca stim from : " , minx, " to " , maxx, " . val = " , CAStim
  for node in ca.nodes:
    if minx < node.x * dend.L < maxx: node.concentration = CAStim
  if ca_stimT > 0: # only re_init if called from event queue
    # rxd.re_init() #
    h.cvode.re_init()

# adjust channel densities at "branch points"
def place_branch (minx=dend.L/2-2, maxx=dend.L/2+2,IP3Origin=20,ERScale=1):
  # print "place_branch: ", minx, maxx, IP3Origin, ERScale
  for node in gip3r.nodes:
    if minx < node.x * dend.L < maxx: node.value = IP3Origin * ERScale
  for node in gleak.nodes:
    if minx < node.x * dend.L < maxx: node.value = node.value * ERScale
  for node in gserca.nodes:
    if minx < node.x * dend.L < maxx: node.value = node.value * ERScale

# set ip3 boosting 
def set_boost (BoostEvery=0):
  if BoostEvery > 0:
    BoostX = dend.L / 2 + BoostEvery
    while BoostX < dend.L:
      place_branch(BoostX - boost_halfw, BoostX + boost_halfw, ip3_origin, er_scale)
      place_branch(dend.L - (BoostX + boost_halfw), dend.L - (BoostX - boost_halfw), ip3_origin, er_scale)
      BoostX += BoostEvery

if ip3_stimT==0: place_ip3_stim(minx=stim_minx,maxx=stim_maxx,IP3Stim=ip3_stim)#called in myrun via event queue otherwise
if ca_stimT==0:  place_ca_stim(minx=stim_minx,maxx=stim_maxx,CAStim=ca_stim)#called in myrun via event queue otherwise

h.dt = dt # only used when variable time-step integration is off (cvode.active()==0)

data = {} # data structure for concentrations - voltArray is for voltage
data["cytca"] = numpy.zeros( (dend.nseg, ceil(simdur/recdt) ) )
data["erca"] = numpy.zeros( (dend.nseg, ceil(simdur/recdt) ) )
data["hgate"] = numpy.zeros( (dend.nseg, ceil(simdur/recdt) ) )
data["ip3"] =  numpy.zeros( (dend.nseg, ceil(simdur/recdt) ) )
datacol = 0

# code for an alternative recording dictionary, using vec.record to reconrd concentrations
def recordCytCa_vec():
  '''use vec.record to record cyt calcium concentration. Will return list containing
  vectors, one vector for each node'''
  cytcaRecordl = []
  for mynode in ca[cyt].nodes:
    cytcavec = h.Vector()
    cytcavec.record(mynode._ref_concentration, recdt)
    cytcaRecordl.append(cytcavec)
  return cytcaRecordl

def recordERCa_vec():
  '''use vec.record to record ER calcium concentration. Will return list containing
  vectors, one vector for each node'''
  ercaRecordl = []
  for mynode in ca[er].nodes:
    ercavec = h.Vector()
    ercavec.record(mynode._ref_concentration, recdt)
    ercaRecordl.append(ercavec)
  return ercaRecordl

def recordHgate_vec():
  '''use vec.record to record h gate "concentration". Will return list containing
  vectors, one vector for each node'''
  hgateRecordl = []
  for mynode in h_gate.nodes:
    hgatevec = h.Vector()
    hgatevec.record(mynode._ref_concentration, recdt)
    hgateRecordl.append(hgatevec)
  return hgateRecordl

def recordIP3_vec():
  '''use vec.record to record ip3 concentration in ER, Will return list containing
  vectors, one vector for each node'''
  ip3Recordl = []
  for mynode in ip3.nodes:
    ip3vec = h.Vector()
    ip3vec.record(mynode._ref_concentration, recdt)
    ip3Recordl.append(ip3vec)
  return ip3Recordl

def recordTI_vec():
  '''use vec.record to record current across T-type of calcium channels. Will return list containing vectors, one vector for each node'''
  TI_Recordl = []
  for mynode in ca[cyt].nodes:
    T_Ivec = h.Vector()
    T_Ivec.record(dend(mynode.x).cat._ref_ica, recdt)
    TI_Recordl.append(T_Ivec)
  return TI_Recordl

def recordNI_vec():
  '''use vec.record to record current across N-type of calcium channels. Will return list containing vectors, one vector for each node'''
  NI_Recordl = []
  for mynode in ca[cyt].nodes:
    N_Ivec = h.Vector()
    N_Ivec.record(dend(mynode.x).can._ref_ica, recdt)
    NI_Recordl.append(N_Ivec)
  return NI_Recordl
  
def recordLI_vec():
  '''use vec.record to record current across N-type of calcium channels. Will return list containing vectors, one vector for each node'''
  LI_Recordl = []
  for mynode in ca[cyt].nodes:
    L_Ivec = h.Vector()
    L_Ivec.record(dend(mynode.x).cal._ref_ica, recdt)
    LI_Recordl.append(L_Ivec)
  return LI_Recordl


data_vec = {}
data_vec['cytca'] = recordCytCa_vec()
data_vec['erca'] = recordERCa_vec()
data_vec['hgate'] = recordHgate_vec()
data_vec['ip3'] = recordIP3_vec()
if electrical:
  data_vec['t_cachan'] = recordTI_vec()
  data_vec['n_cachan'] = recordNI_vec()
  data_vec['l_cachan'] = recordLI_vec()

timeVec = h.Vector(); timeVec.record(h._ref_t)

# record voltage using a reference to voltage ._ref_v using rxd species node locations
def recordVolt_vec():
  '''record dendritic membrante voltage'''
  voltlist = []
  for node in ca[cyt].nodes(dend):
    myvec = h.Vector()
    myvec.record(dend(node.x)._ref_v, recdt)
    voltlist.append(myvec)
  return voltlist

voltlist = recordVolt_vec()

# save data dictionary as compressed numpy array (npz format)
def mysavedata (simstr,ldata=data):
  fname = "./data/" + simstr + "_.npz"
  print 'saving data...'
  if electrical:
    numpy.savez_compressed(fname,cytca=ldata["cytca"],erca=ldata["erca"],hgate=ldata["hgate"],ip3=ldata["ip3"], volt=np.array(voltlist))
  else:
    numpy.savez_compressed(fname,cytca=ldata["cytca"],erca=ldata["erca"],hgate=ldata["hgate"],ip3=ldata["ip3"])
  
# load data dictionary
def myloaddata (simstr):
  fname = "./data/" + simstr + "_.npz"
  ldata = {}
  xn = numpy.load(fname)
  for k in xn.keys(): ldata[k] = xn[k]
  del xn.f
  xn.close()
  return ldata

# initialize recording data structures
def initrec ():
  global datacol
  datacol = 0
  for k in data.keys(): data[k] = numpy.zeros( (dend.nseg, ceil(simdur/recdt))  )

# record Nodelist nl's concentrations into dat (2D numpy array) 
def dorec (dat, nl):
  i = 0 #to go over each node
  for node in nl:
    dat[i][datacol] = node.concentration
    i += 1

# record cytoplasmic ca concentration
def recCYTca (): dorec(data["cytca"], ca[cyt].nodes(dend))

# record ER ca concentration
def recERca (): dorec(data["erca"], ca[er].nodes(dend))

# record IP3 hgate values
def recHGate (): dorec(data["hgate"], h_gate.nodes(dend))

# record IP3 concentration
def recIP3 (): dorec(data["ip3"], ip3.nodes(dend))

IL_loclist = [] # a list of vectors used to record current accross L-channels
gL_loclist = [] # a list to record conductance accross L-channels
ampaCurrl = [] # list for AMPA current
if electrical:
  for node in ca[cyt].nodes:
    myvec = h.Vector()
    myvec.record(dend(node.x).cal._ref_ica, recdt)
    IL_loclist.append(myvec)
  for node in ca[cyt].nodes:
    myvec = h.Vector()
    myvec.record(dend(node.x).cal._ref_gcal, recdt)
    gL_loclist.append(myvec)
  for mysyn in synlist:
    myvec = h.Vector()
    myvec.record(mysyn._ref_i, recdt)
    ampaCurrl.append(myvec)

# records calcium and voltage- called during run
def dorecord ():
  global datacol
  recCYTca()
  recERca()
  recHGate()
  recIP3()
  datacol += 1

# run the sim and save calcium levels in data
def myrun ():
  clockStart = clock()
  if h.t > 0: # was this sim run already?
    if ip3_stimT==0: place_ip3_stim(minx=stim_minx,maxx=stim_maxx,IP3Stim=ip3_stim)
    if ca_stimT==0: place_ca_stim(minx=stim_minx,maxx=stim_maxx,CAStim=ca_stim)  
  print "starting simulation..."
  initrec() # initialize recording data structures - does not create any events or call Vector record
  h.init() # contains a call to h.finitialize
  # must setup hotspots after h.init/h.finitialize, otherwise the rxd.Parameter initial values are set/used
  place_branch(minx=dend.L/2-boost_halfw,maxx=dend.L/2+boost_halfw,IP3Origin=ip3_origin,ERScale=er_scale)
  set_boost(BoostEvery=boost_every)
  if loadState:
    loadstate(statestr) # load initial values from saved file (hoc/rxd)
  elif useInitDict:
    setHocStateVars(fcfg) # load initial values from config file
    setRxDStateVars(fcfg)
  else: # make sure ca[er].concentration initialized
    cae_init = (caAvg_init - caCYT_init * fc) / fe
    ca[er].concentration = cae_init
    h.cvode.re_init()
  if IP3ForceInit:
    ip3[cyt].concentration = ip3_init
    h.cvode.re_init() 
  print '\ntstart:',tstart,' tstop:', tstop, 'h.tstop:', h.tstop
  # only set event if supposed to place stim(s) after t=0 (otherwise already set above)
  if ip3_stimT > 0: h.cvode.event(tstart+ip3_stimT,place_ip3_stim) # put an event on queue to set ip3 stim
  if ca_stimT > 0: h.cvode.event(tstart+ca_stimT,place_ca_stim) # put an event on queue to set ca stim
  for t in numpy.arange(tstart, tstop, recdt): h.cvode.event(t, dorecord)
  for t in numpy.arange(tstart, tstop, 500): h.cvode.event(t, displaySimTime)#displays h.t every 500 msec
  # set synapse properties based on specification in config file - only changing number,active seems to work
  # using FInitializeHandler type 2,3 did not help
  if electrical: setSynProps() 
  h.continuerun(tstop) # does not guarantee that tstop is reached
  while h.t < tstop: h.advance() # make sure get to tstop
  #print 'ran from ' , timeVec[0], ' to ', timeVec[1], ' to ', timeVec[-1]  
  #printt('final time')  
  clockEnd = clock()  
  print '\nsim runtime:',str(datetime.timedelta(seconds=clockEnd - clockStart)),'secs'#format s as hh:mm:ss

# loop using fadvance and save the dts - only useful when have cvode.active(1)
def checkdts (iters=1000):
  ldt,lt = [],[]
  for i in xrange(iters):
    ldt.append(h.dt)
    lt.append(h.t)
    h.fadvance()
  return ldt,lt

# draw calcium level using dat (2D numpy array)
def mydraw (data=data,vmin=[0,0,0],vmax=[0.002,0.2, 0.011],startt=tstart,endt=tstop,\
            miny=0,maxy=int(dend.L),keys=["cytca"],titles=["ca[cyt] (mM)","ip3[cyt] (mM)"],\
            interp='bilinear',mycmap=cm.jet):
  #ion() # turn on interactive
  mp.figure()
  gn,i = 1, 0
  splotlist = []
  for k in keys:
    dat = data[k]
    splot = subplot(len(keys),1,gn);  mymin,mymax = amin(dat),amax(dat);
    imshow(dat, vmin=vmin[i], vmax=vmax[i],aspect='auto',extent=(startt/1e3,endt/1e3,0,int(dend.L)), origin='lower',interpolation=interp,cmap=mycmap)
    xlim( (startt/1e3, endt/1e3) ); ylim( (miny, maxy+1) );
    ylabel(r'Position ($\mu$m)');  title(titles[i]); colorbar();
    if i == len(keys)-1: xlabel('Time(s)');
    print keys[i],":",mymin, mymax
    splotlist.append(splot)
    gn += 1
    i += 1
  return splotlist


def markStim(mysubplot):
  '''will mark time and location of stimulation on subplots generated by mydraw'''
  ip3stim_x = np.array(get_ip3stimLoc())*dend.L
  ip3stim_time = np.repeat(ip3_stimT/1e3, len(ip3stim_x))
  mysubplot.plot(ip3stim_time, ip3stim_x, 'yo', markersize=8)
  mysubplot.title.set_text(mysubplot.title.get_text()+' [circle: IP3 stim]')#, cross: syn depol.]')
  if len(netConRecveclist)>0:
    mysubplot.title.set_text(mysubplot.title.get_text()+' [cross: syn depol.]')
    for myvec in enumerate(netConRecveclist):
      #debug_here()
      Lchan_stim_time = np.array(myvec[1])/1e3#; print Lchan_stim_time
      Lchan_stim_x = np.repeat(synLocl[myvec[0]]*dend.L, Lchan_stim_time.size)#; print Lchan_stim_x
      mysubplot.plot(Lchan_stim_time, Lchan_stim_x, 'mx', markersize=12, markeredgewidth=2)


def voltDraw(voltArr=np.array(voltlist), startt=tstart, endt=h.tstop, miny=0, maxy=int(dend.L), mycmap=cm.jet):
  '''will plot voltage from voltArray (2D image)'''
  mp.gca()
  debug_here()
  mymin,mymax = amin(voltArr),amax(voltArr)
  mp.imshow(voltArr, aspect='auto', extent=(startt/1e3, endt/1e3, 0, int(dend.L)), cmap=mycmap)
  ylabel(r'Position ($\mu$m)');   xlabel('Time (s)'); title('Voltage (mV)')
  mp.colorbar()


def testDefaultParams(voltArr=np.array(voltlist)):
  '''test default params'''
  return voltArr

# apply interpolation to x,y (kind specifies method, cubic, linear, etc.)
def csmoothxy (x,y,sz=0,kind='cubic'):
  fs = interp1d(x,y,kind=kind)
  if sz < 1: sz = 3*len(x)
  xnew = numpy.linspace(amin(x), amax(x), sz)
  ynew = fs(xnew)
  return xnew,ynew

# cut out the individual waves via thresholding and component labeling
def wavecut (im,thresh=caCYT_init*2):
  mask = im > thresh
  labelim, nlabels = ndimage.label(mask)
  return labelim, nlabels

# scan horizontally
def scanx (lab,thresh,widx,x,y,slicex,yh,lastyh):
  endx,speed = x,0
  while endx < slicex.stop and lab[y][endx] == widx: endx += 1
  if endx >= slicex.stop: endx = slicex.stop - 1
  # yh and lastyh are in units of pixels. to convert to um/ms need to use spaceum (um / segment)
  if lastyh >= 0: speed = spaceum*(yh - lastyh) / recdt # um / ms
  return endx,speed

# get the wave information in an nqs
#  extraction based on threshold crossing, then image labeling (connected comps)
#  NB: speed calculation assumes that IP3 starts at bottom and moves upwards.
#  output nqs columns: widx=wave id, startt=start time (ms), endt=end time (ms)
#   y=position on dend, speed=speed from last position in um/ms, durt=duration (ms)
#   overlap=whether run overlaps with prior run (run=horizontal segment of a wave)
#   up=whether was looking up (1) or down (0) for the particular run
def wavenq (dat,thresh=caCYT_init*2,verbose=False):
  nq = NQS("widx","startt","endt","y","speed","durt","overlap","up")
  lab, nwaves = wavecut(dat,thresh)
  nq.clear(len(dat[0]*nwaves))
  for widx in xrange(1,nwaves+1):
    slicey, slicex = ndimage.find_objects(lab==widx)[0]
    midy = int(slicey.start + (slicey.stop - slicey.start) / 2.0)
    lastyup,lastydown,speed = midy,midy,0
    lastyuph,lastydownh = midy,midy
    if verbose:
      print "slicex:",slicex.start, slicex.stop
      print "slicey:",slicey.start, slicey.stop   
    # y0 + (y1-y0) * (threshold-f(y0)) /  (f(y1)-f(y0))
    x,y = slicex.start, midy
    if lab[y][x] == widx:
      endx,speed = scanx(lab,thresh,widx,x,y,slicex,y,lastyup)
      nq.append(widx,x*recdt,endx*recdt,y*spaceum,speed,recdt*(endx-x),0,1)
    lastendxup,lastendxdown = -1,-1
    while x < slicex.stop: # traverse through time      
      found = False 
      y = slicey.stop - 1 # look for highest point
      if verbose: print "x:",x, "y:", y, " = ", lab[y][x]
      while y >= slicey.start: # starting above and going down until hit it
        if lab[y][x] == widx:
          if x > 0:
            if lab[y][x-1] != widx: # make sure on outer edge (avoids overlap)
              found = True
              break
          else:
            found = True
            break
        y -= 1        
      if found:
        yh = y
        if y + 1 < len(dat):
          y0,y1 = y, y + 1 # y0 is part of wave, y1 is above its top at t=x so dat[y1][x] <= thresh
          yh = y0 + (thresh-dat[y1][x]) / (dat[y0][x]-dat[y1][x])
        if verbose: print "found u " , x , y, yh
        endx,speed = scanx(lab,thresh,widx,x,y,slicex,yh,lastyuph)
        olap = 0
        if y == lastyup and x < lastendxup: olap = 1
        nq.append(widx,x*recdt,endx*recdt,yh*spaceum,speed,recdt*(endx-x),olap,1);
        lastyup, lastendxup, lastyuph = y, endx, yh
      found = False 
      y = slicey.start # look for lowest point
      while y < slicey.stop: # starting below and going up until hit it
        if lab[y][x] == widx:
          if x > 0:
            if lab[y][x-1] != widx: # make sure on outer edge (avoids overlap)
              found = True
              break
          else:
            found = True
            break
        y += 1
      if found:
        yh = y
        if y - 1 >= 0: 
          y0,y1 = y, y - 1 # y0 is part of wave, y1 is below its bottom at t=x so dat[y1][x] <= thresh
          yh = y0 - (thresh-dat[y1][x]) / (dat[y0][x]-dat[y1][x])
        if verbose: print "found d " , x , y, yh
        endx,speed = scanx(lab,thresh,widx,x,y,slicex,yh,lastydownh)
        olap = 0
        if y == lastydown and x < lastendxdown: olap = 1
        nq.append(widx,x*recdt,endx*recdt,yh*spaceum,speed,recdt*(endx-x),olap,0);
        lastydown, lastendxdown, lastydownh = y, endx, yh
      x += 1 # move to next time-point
  return nq

# draw cytca wave with overlay of start/end/top
def mydrawwave (data=data,key="cytca",nqw=None,first=True,last=True):
  dodel = False
  if nqw is None:
    nqw = wavenq(data[key])
    dodel = True
  mydraw(data,keys=[key])
  dw = nqs2pyd(nqw)
  if first: plot(dw["startt"]/1e3,dw["y"],'ko',linewidth=3)
  if last: plot(dw["endt"]/1e3,dw["y"],'ro',linewidth=3)
  if dodel: nqsdel(nqw)

# get first wave to occur after the stim time (stimt). nqw is from getwavenq
def firstwaveaftert (nqw,stimt):
  nqw.verbose=0; nqw.tog("DB")
  mxw,idx=int(nqw.getcol("widx").max()),-1
  for widx in xrange(1,mxw+1,1):
    nqw.select("widx",widx)
    st=nqw.getcol("startt").min()
    #print widx, st
    if st >= stimt:
      idx = widx
      break
  nqw.verbose=1
  return idx

# make a dictionary from nqw
def makedw (data,key,nqw,nooverlap,NM=1,firstaftstim=True):
  dodel = False
  if nqw is None:
    nqw = wavenq(data[key])
    dodel = True
  widx = 1
  if firstaftstim: widx = firstwaveaftert(nqw,ip3_stimT)
  nqw.verbose,N=0,0
  if nooverlap:
    N = nqw.select("y",">=",500,"up",1,"widx",widx,"overlap",0) # only want upward portion of 1st wave
  else:
    N = nqw.select("y",">=",500,"up",1,"widx",widx) # only want upward portion of 1st wave
  debug_here()# value of N, NM???
  if N < NM:
    nqw.verbose=1
    if dodel: nqsdel(nqw)
    return 0.0
  dw = nqs2pyd(nqw)
  nqw.verbose=1
  if dodel: nqsdel(nqw)
  return dw

# get the wave speed in um/s (instantaneous estimate)
def instwavespeed (data=data,key="cytca",nqw=None,nooverlap=True):
  dw = makedw(data,key,nqw,nooverlap,NM=2)
  if dw == 0: return 0.0
  sp = [0]
  for i in xrange(1,len(dw["y"]),1):
    dy = dw["y"][i]-dw["y"][i-1]
    dx = dw["startt"][i]-dw["startt"][i-1]
    if dx > 0:
      sp.append( 1e3*spaceum*dy/dx )
    else:
      sp.append( 0 )
  return numpy.array(sp)

# get the wave speed in um/s (uses first and last position)
def wavespeed (data=data,key="cytca",nqw=None,nooverlap=True):
  #debug_here()
  dw = makedw(data,key,nqw,nooverlap,NM=2)
  if dw == 0: return 0.0
  try:
    dy = dw['y'][-1]-dw['y'][0]
    dx = dw['startt'][-1]-dw['startt'][0]
    if dx > 0.0:
      return 1e3 * dy / dx
    else:
      return 0.0
  except:
    return 0.0

# get the wave duration
def wavedur (data=data,key="cytca",nqw=None,func=numpy.median,nooverlap=True):
  dw = makedw(data,key,nqw,nooverlap,NM=1)
  if dw == 0.0: return 0.0
  return func( dw["durt"] )

# get the full duration of the 1st wave (last time - first time)
def wavetime (data=data,key="cytca",nqw=None,nooverlap=True):
  dw = makedw(data,key,nqw,nooverlap,NM=1)
  if dw == 0.0: return 0.0
  return amax( dw["endt"] ) - amin( dw["startt"] )

# get the max position wave has traveled 
def wavedist (data=data,key="cytca",nqw=None,nooverlap=True):
  dw = makedw(data,key,nqw,nooverlap,NM=1)
  if dw == 0.0: return 0.0
  return amax( dw['y'] )

# get the wave's starting time
def waveonset (data=data,key="cytca",nqw=None,nooverlap=True):
  dw = makedw(data,key,nqw,nooverlap,NM=1)
  if dw == 0.0: return 0.0
  mint = amin(dw["startt"])
  return mint

# save figure
def mysavefig ():
  hg_version = subprocess.check_output(['hg', 'identify', '-i']).strip()
  outty = ".png"
  outf = "gif/cawave_{0}_{1}_{2}_{3}_{4}_{5}_{6}" + outty
  savefig(outf.format(hg_version, er_scale, ip3_notorigin, ip3_origin, gserca0, boost_every, ip3_stim), transparent=True)

# run and draw
def myrundraw ():
  '''run and draw'''
  myrun()
  mydraw()

# run, draw, save! 
def myrundrawsv ():
  myrundraw()
  mysavefig()

# saves NMODL/non-rxd state variables
def savehocstate(filestr):
  '''saves the current hoc states (not rxd) of the simulation into file with name filestr.dat
  filestr: a string of file path (without file extension)'''
  myss = h.SaveState()
  #h.initnrn() # to set h.t to zero
  myss.save()
  myfile = h.File('./data/'+filestr+'_hoc_.dat')
  myss.fwrite(myfile)
  myfile.close()

# restores NMODL/non-rxd state variables and h.t  
def loadhocstate (filestr):
  '''loads hoc states (not rxd states) of simulation from file filestr.dat
  filestr: a string of file path (without file extension).
  For loading to happen appropriately, this statement has to be called after any initialization (eg h.finitalize)'''
  global tstart, tstop
  myfile = h.File()
  myfile.ropen('./data/'+filestr+'_hoc_.dat')
  myss = h.SaveState()
  if myfile.isopen():
    myss.fread(myfile)
    myfile.close()
    myss.restore(1)
    print 'loaded hoc states from:', myfile.getname()
    tstart = h.t
    tstop = tstart + simdur # stopping time of simulation
    h.tstop = tstop
    #for ns in netStimlist: ns.start = nstimStart + tstart
  else: print "file cannot be open to read"

# save rxd-related state variables
def saverxdstate (filestr):
  '''save rxd state into file filestr.npy'''
  try:
    np.save('./data/'+filestr+'_rxd_.npy', rxd.node._states)
  except:
    print 'no rxd'

# load rxd-related state variables
def loadrxdstate(filestr):
  '''load rxd state
  for loading to occur appropriately, it has to take place after any initialization stattment (eg. h.finitialize)'''
  fname = './data/'+filestr+'_rxd_.npy'
  xx = None
  try:
    xx = np.load(fname)
  except:
    print 'loadrxdstate ERRA: could not load ', fname
  try:
    rxd.node._states[:] = xx # this line has to be that way to avoid memory problems
    print 'loaded rxd states from:', fname
  except:
    print 'loadrxdstate ERRB: could not set rxd.node._states from ', fname

# saves simulation state and time
def savestate (filestr):
  '''saves both hoc and rxd states into files with extensions filestr.dat & filestr.npy respectively'''
  savehocstate(filestr)
  saverxdstate(filestr)

def loadstate (filestr):
  '''loads hoc and rxd states from files with extensions filestr.dat & filestr.npy, respectively'''
  if not loadRXDStateOnly: loadhocstate(filestr)
  loadrxdstate(filestr)
  h.cvode.re_init()

def displaySimTime():
  '''displays h.t as the simulation is running, as numbers that are changing dynamically - helpful to keep track of simulations that are running for long'''
  sys.stdout.write('\rh.t: {0} msec...'.format(h.t))
  sys.stdout.flush()

#guibox = makegui() # delete this if not needed 

if __name__ == '__main__':     # if ran directly
  if runit: # run sim ?
    print "running..."
    myrun()
    if saveout: # save data ?
      print "saving output..."
      mysavedata(simstr,ldata=data)
    if saveState: # save state info
      print 'saving state...'
      savestate(simstr) # save the state to simstr+.dat (for hoc) and simstr+.npy (for rxd)
    if dodraw:
      print "drawing..."
      mydraw()
      show()
