from neuron import h
h.load_file("nqs.hoc")
import numpy
from vector import *

NQS = h.NQS
nqsdel = h.nqsdel

# converts 2D numpy array into an NQS
def np2nqs (npa,names=[]):
  nrow,ncol = numpy.shape(npa)
  if ncol < 2:
    print "np2nqs ERRA: must have at least 2 columns!"
    return None
  nqo = NQS(ncol)
  for i in xrange(ncol):
    vcol = py2vec( npa[0:nrow][i] )
    nqo.v[i].copy(vcol)
  for i in xrange(len(names)): # assign names
    nqo.s[i].s = names[i]
  return nqo

# converts nqs to numpy array
def nqs2np (nq,Full=False):
  if Full: nq.tog("DB")
  ncol,nrow = int(nq.m[0]),nq.size()
  npa = numpy.zeros( (nrow,ncol) )
  for i in xrange( ncol ):
    nptmp = vec2nparr( nq.getcol(nq.s[i].s) )
    npa[0:nrow][i] = nptmp[0:nrow]
  return npa

# convert nqs to a python dictionary
def nqs2pyd (nq,Full=False):
  d = {};
  if Full: nq.tog("DB");
  for i in xrange(int(nq.m[0])): d[nq.s[i].s] = vec2np(nq.getcol(nq.s[i].s))
  return d

