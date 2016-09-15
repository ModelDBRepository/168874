NEURON { SUFFIX naf }
NEURON {  USEION na WRITE ina }
ASSIGNED { ina }
PARAMETER {
	erev 		= 60.  (mV)
	gmax 		= 0.030    (mho/cm2)
        vrest           = -60.

	exptemp		= 37
	maflag 		= 3
	malphaA 	= -0.32
	malphaB		= -4.0
	malphaV0	= 13.1
	mbflag 		= 3
	mbetaA 		= 0.28
	mbetaB		= 5.0
	mbetaV0		= 40.1
	mq10		= 3
	mexp 		= 2

	haflag 		= 1
	halphaA 	= 0.128
	halphaB		= -18
	halphaV0	= 22.
	hbflag 		= 2
	hbetaA 		= 4.
	hbetaB		= -5.
	hbetaV0		= 45.
	hq10		= 3
	hexp 		= 1

	vmax 		= 100  (mV)
	vmin 		= -100 (mV)
} : end PARAMETER

INCLUDE "geneval_cvode.inc"

PROCEDURE iassign () { i = g*(v-erev) ina=i }
