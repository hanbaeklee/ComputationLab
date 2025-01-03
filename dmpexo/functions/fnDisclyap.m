function V = fnDisclyap(a1,b1,vecflag)
% disclyap -- extension of Hansen/Sargent's DOUBLEJ.M
%
%  function V = disclyap(a1,b1,vecflag)
%  Computes infinite sum V given by
%
%         V = SUM (a1^j)*b1*(a1^j)'
%
%  where a1 and b1 are each (n X n) matrices with eigenvalues whose moduli are
%  bounded by unity, and b1 is an (n X n) matrix.
%  The sum goes from j = 0 to j = infinity.  V is computed by using
%  the following "doubling algorithm".  We iterate to convergence on
%  V(j) on the following recursions for j = 1, 2, ..., starting from
%  V(0) = b1:
%
%       a1(j) = a1(j-1)*a1(j-1)
%       V(j) = V(j-1) + a1(j-1)*V(j-1)*a1(j-1)'
%
%  The limiting value is returned in V.
%
% -----------------------------------------------------------------------------
% EMT added following comments and the vecflag argument (default=false)
% -----------------------------------------------------------------------------
%
% 1) L&S p. 76, write
% X = doublej(A, B * B') solves X = A * X * A' + B * B' where a1 = A and B*B' = b1
%
% 2) Note on algorithm
% * let $V_k \equiv \sum_{j=0}^K A^j C'C A^j'$
% * then $V_{2k} = V_k + A^k V_k A^k'$
% * let $k=2^j$ and rewrite $V_k \equiv V_j$ and $A_j\equiv A^{2^j}$ so that $A_{j+1} = A_j^2$
% * then $V_{j+1} = V_j + A_j^2 V_j A_j^2'$
%
% the doubling uses fewer iterations (namely only $j$ instead of $2^j$ and is thus much faster
%
% 3) Solves "discrete Lyapunov" equation (but can also handle constants in the SS)
%
% 4) new parameter: if vecflag == true: apply vec formula of Hamilton instead of doubling
%    will assume invertibility, i.e. no constant in SS

error(nargchk(2,3,nargin))
if nargin < 3
   vecflag = false;
end

if vecflag
   Np          = size(a1, 1);
   Inp2        = eye(Np^2);
   vecSIGMA    = ((Inp2 - kron(a1, a1)) \ Inp2) * b1(:);
   V           = reshape(vecSIGMA, Np, Np);
else
   alpha0   =  a1;

   gamma0   =  b1;

   delta  =  5;
   ijk   =  1;

   while delta > eps;

      alpha1   =  alpha0 * alpha0;
      gamma1   =  gamma0 + alpha0 * gamma0 * alpha0';
      %    delta    =  max(max(abs(gamma1-gamma0)));
      delta    = max(abs(gamma1(:) - gamma0(:)));
      gamma0   = gamma1;
      alpha0   = alpha1;

      ijk      = ijk + 1;

      if ijk > 50;
         error('Error: ijk = %d, delta = %f check aopt and c for proper configuration', ijk, delta)
      end

   end

   V = gamma1;
end