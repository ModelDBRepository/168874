import ConfigParser
import io
import sys

# default config as string
def_config = """
[set]
boost_every = 0.0
boost_halfw = 5.0
er_scale = 1.0
ca_stim = 0.0
ca_stimT = 2000.0
ip3rtau = 400
ip3_notorigin = 120400.0
ip3_origin = 120400.0
gserca = 1.9565
gleak = 18.06
stim_minx = 498
stim_maxx = 502
ip3_stim = 1.25
ip3_stimT = 2000
IP3ForceInit = 0
ip3_init = 0.1
cacyt_init = 0.0001
caAvg_init = 0.0017
cadiff = 0.080
ip3diff = 1.415
gCaChannels = 1.e-6
spaceum = 1.0
synLoclist = [0.5]
nstimStart = 1000
nstimInterval = 1
nstimNumber = 40
nconnThreshold = 0
nconnDelay = 0
nconnWeight = 20
nconnActive = 1
electrical = 0
[run]
recdt = 5
tstop = 30000
dt = 1
cvodeactive = 1
runit = 0
simstr = '14aug5_A'
saveout = 0
dodraw = 0
loadState = 0
ipydebug = 0
statestr = 'mystate'
savestate = 0
useInitDict = 0
loadRXDStateOnly = 1
"""

# write config file starting with defaults and new entries
# specified in section (sec) , option (opt), and value (val)
# saves to output filepath fn
def writeconf (fn,sec,opt,val):
  conf = ConfigParser.ConfigParser()
  conf.readfp(io.BytesIO(def_config)) # start with defaults
  # then change entries by user-specs
  for i in xrange(len(sec)): conf.set(sec[i],opt[i],val[i])
  # write config file
  with open(fn, 'wb') as cfile: conf.write(cfile)

# read config file
def readconf (fn="cawave.cfg"):

  config = ConfigParser.ConfigParser()
  config.read(fn)

  def conffloat (base,var,defa): # defa is default value
    val = defa
    try: val=config.getfloat(base,var)
    except: pass
    return val

  def confint (base,var,defa):
    val = defa
    try: val=config.getint(base,var)
    except: pass
    return val

  def confstr (base,var,defa):
    val = defa
    try: val = config.get(base,var)
    except: pass
    return val
  
  def conflist (base, var, defa):
    val = defa
    try: val = config.get(base, var)
    except: pass
    return eval(val) # eval will convert the string to a list

  er_scale = conffloat("set","er_scale", 1.0)
  ip3_notorigin = conffloat("set","ip3_notorigin", 120400.0)
  ip3_origin = conffloat("set","ip3_origin", 120400.0)
  gserca0 = conffloat("set","gserca", 1.9565)
  boost_every = conffloat("set","boost_every", 0.0)
  boost_halfw = conffloat("set","boost_halfw", 5.0)
  gleak0 = conffloat("set","gleak", 18.06)
  ip3rtau = conffloat("set","ip3rtau", 400.0)
  recdt = conffloat("run","recdt", 5.0)
  tstop = conffloat("run","tstop", 30000.0)
  cvodeactive = confint("run","cvodeactive", 1)
  dt = conffloat("run","dt", 50.0)
  runit = conffloat("run","runit", 0.0)
  simstr = confstr("run","simstr","13nov18_")
  saveout = conffloat("run","saveout",0.0)
  stim_minx = conffloat("set","stim_minx",123)#498.0)
  stim_maxx = conffloat("set","stim_maxx",127)#502.0)
  dodraw = conffloat("run","dodraw",0.0)
  loadState = confint('run', 'loadState', 0)
  loadRXDStateOnly = confint('run','loadRXDStateOnly',1)
  ipydebug = confint('run','ipydebug',0)
  savestate = confint("run","savestate",0)
  useInitDict = confint("run","useInitDict",0) # whether to use sstvars init dict
  statestr = confstr('run', 'statestr', 'mystate')
  ip3_stim = conffloat("set","ip3_stim", 1.25)
  ip3_stimT = conffloat("set","ip3_stimT",2000.0)
  IP3ForceInit = confint("set","IP3ForceInit",0)
  ip3_init = conffloat("set","ip3_init",0.1)
  caCYT_init = conffloat("set","caCYT_init",0.0001)
  caAvg_init = conffloat("set","caAvg_init",0.0017)
  caDiff = conffloat("set","caDiff",0.080) # calcium diffusion coefficient
  ip3Diff = conffloat("set","ip3Diff",1.415) # ip3 diffusion coefficient
  gCaChannels = conffloat('set', 'gCaChannels', 1.e-6)
  spaceum = conffloat("set","spaceum",1.0) # spatial resolution for each segment, in um (nseg=L/spaceum)
  ca_stim = conffloat("set","ca_stim",0.0)
  ca_stimT = conffloat("set","ca_stimT",2000.0)
  synLoclist = conflist('set', 'synLoclist', []) # list of locations to insert depolarizing synapses
  nstimStart = conffloat('set', 'nstimStart', 1000)
  nstimInterval = conffloat('set', 'nstimInterval', 1)
  nstimNumber = conffloat('set', 'nstimNumber', 40)
  nconnThreshold = conffloat('set', 'nconnThreshold', 0)
  nconnDelay = conffloat('set', 'nconnDelay', 0)
  nconnWeight = conffloat('set', 'nconnWeight', 20)
  nconnActive = confint('set', 'nconnActive', 1)
  electrical = confint('set', 'electrical', 0) # whether to include ion channels and synapses
  
  sys.stdout = open("./simconfloaded.log", "w") # will direct the print output to simconfloaded.log file
  print "er_scale:",er_scale,",ip3_notorigin:",ip3_notorigin,",ip3_origin:",ip3_origin,\
      ",gserca:",gserca0,",boost_every:",boost_every,",ip3_stim:",ip3_stim, ",ip3rtau:",ip3rtau,\
      "cvodeactive:",cvodeactive," ,dt:",dt," ,runit:",runit," ,simstr:",simstr," ,saveout:",saveout,\
      "stim_minx:",stim_minx," ,stim_maxx:",stim_maxx," ,dodraw:",dodraw, ', loadState:',loadState, ", ip3_stimT:",ip3_stimT,\
      " ,ip3_init:",ip3_init, " ,caCYT_init:",caCYT_init, " ,caAvg_init:",caAvg_init,\
      " ,caDiff:",caDiff, " ,ip3Diff:",ip3Diff, " ,gCaChannels:",gCaChannels, " ,spaceum:",spaceum,\
      " ,ca_stim:",ca_stim, " ,ca_stimT:",ca_stimT," ,recdt:",recdt,\
      " ,boost_halfw:",boost_halfw," ,gleak0:",gleak0, ' ,synLoclist:', synLoclist, ' ,nstimStart:', nstimStart,\
      ' ,nstimInterval:', nstimInterval, ' ,nstimNumber:', nstimNumber, ' ,nconnThreshold:', nconnThreshold,\
      ' ,nconnDelay:', nconnDelay, ' ,nconnWeight:', nconnWeight, ' ,nconnActive:', nconnActive, \
      ' ,ipydebug:', ipydebug,  ' ,statestr:', statestr, ' ,savestate:', savestate, ' ,useInitDict:', useInitDict, \
      ' ,IP3ForceInit:', IP3ForceInit, ' ,loadRXDStateOnly:', loadRXDStateOnly, ' ,electrical:', electrical
  sys.stdout = sys.__stdout__
  return er_scale,ip3_notorigin,ip3_origin,gserca0,boost_every,ip3_stim,\
      gleak0,ip3rtau,recdt,tstop,cvodeactive,dt,runit,simstr,saveout,\
      stim_minx,stim_maxx,dodraw,loadState,ip3_stimT,ip3_init,caCYT_init,caAvg_init,\
      caDiff,ip3Diff,gCaChannels,spaceum,ca_stim,ca_stimT,boost_halfw,synLoclist,\
      nstimStart, nstimInterval, nstimNumber, nconnThreshold, nconnDelay, nconnWeight, nconnActive, \
      ipydebug, statestr, savestate, useInitDict, IP3ForceInit, loadRXDStateOnly, electrical
