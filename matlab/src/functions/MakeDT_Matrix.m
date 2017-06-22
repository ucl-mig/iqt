function dt = MakeDT_Matrix(d11, d12, d13, d22, d23, d33)
% Makes a diffusion tensor matrix out of the six elements.
%
% dt = MakeDT_Matrix(d11, d12, d13, d22, d23, d33)
% returns the matrix
% [d11 d12 d13]
% [d12 d22 d23]
% [d13 d23 d33]
%
% $Id: MakeDT.m,v 1.1 2007/07/25 09:54:19 ucacdxa Exp $
% 
% ---------------------------
% Part of the IQT matlab package
% https://github.com/ucl-mig/iqt
% (c) MIG, CMIC, UCL, 2017
% License: LICENSE
% ---------------------------
%

dt = [d11 d12 d13; d12 d22 d23; d13 d23 d33];

