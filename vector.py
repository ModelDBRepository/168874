from neuron import h
import numpy

Vector = h.Vector

# turn a NEURON Vector into a numpy array
def vec2np (vec):
  pl = vec.to_python()
  pl = numpy.array(pl)
  return pl

# return Vector from python list or numpy array
def py2vec (pl):
  vec = Vector()
  vec.from_python(pl)
  return vec

# create a random , uniform vector
def RandUnifVec (rdmS,sz,amp):
  vrd = Vector(sz)
  rdm = h.Random()
  rdm.ACG(rdmS)
  rdm.uniform(-amp,amp)
  vrd = Vector(sz)
  vrd.setrand(rdm)
  return vrd

# create a random, normal vector. rdmS=seed, sz=size of vector, mu=average, std=standard-dev
def RandNormVec (rdmS,sz,mu,std):
  vrd = Vector(sz)
  rdm = h.Random()
  rdm.ACG(rdmS)
  rdm.normal(mu,std**2) # takes mean and variance (std**2)
  vrd = Vector(sz)
  vrd.setrand(rdm)
  return vrd

