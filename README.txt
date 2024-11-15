------------------------------------------------------------------

Demo software for Compressive sensing via L1-L2

Reference:

1. Yifei Lou, Penghang Yin, Qi He and Jack Xin, Computing Sparse Representation in a Highly Coherent Dictionary Based on Difference of L1 and L2. to appear in J. of Sci. Computing, 2014. (published online Oct 2014)

2. Penghang Yin, Yifei Lou, Qi He and Jack Xin, Minimization of L1-L2 for Compressed Sensing. to appear in SIAM J. of Sci. Computing, 2015.


------------------------------------------------------------------

Copyright (c) Yifei Lou 
https://sites.google.com/site/louyifei/
This work should be used only for nonprofit purposes.


------------------------------------------------------------------
Contents
------------------------------------------------------------------

The package comprises these functions

*) constrainedL1.m 		: constrained L1 formulation via ADMM
*) constrainedL1L2.m		: constrained L1-L2 formulation via DCA+ADMM
*) unconstrainedL1L2.m		: unconstrained L1-L2 formulation via DCA+ADMM
*) lasso.m			: unconstrained L1 formulation via ADMM 
   (thanks to http://web.stanford.edu/~boyd/papers/admm/lasso/lasso.html)

*) coherence.m			: compute the coherence of a matrix
*) randsample_separated.m	: generate a random sparse vector with separated spikes

*) demo.m			: a demo code 



------------------------------------------------------------------
Feedback
------------------------------------------------------------------

If you have any comment, suggestion, or question, please do
contact Yifei Lou at: yflou@unc.edu

