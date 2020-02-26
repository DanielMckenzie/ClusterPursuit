function Sest = subspacepursuit(Phi,u,K,tol,maxiterations)
%        Sest = subspacepursuit(Phi,u,K,tol,maxiterations)
% Subspace Pursuit algorithm
%   Input
%       K : sparsity of Sest
%       Phi : measurement matrix
%       u: measured vector
%       tol : tolerance for approximation between successive solutions. 
%   Output
%       Sest: Solution found by the algorithm
%
% Algorithm as described in "CoSaMP: Iterative signal recovery from 
% incomplete and inaccurate samples" by Deanna Needell and Joel Tropp.
% 
% This implementation was adapted from a version of CoSaMP by David Mary, 
% but modified and corrected 20110707 by Bob L. Sturm
% Modified again by Daniel Mckenzie June/July 2017 to find least square
% solution approximately using MATLAB's lsqr, instead of exactly.
% 
%
% This script/program is released under the Commons Creative Licence
% with Attribution Non-commercial Share Alike (by-nc-sa)
% http://creativecommons.org/licenses/by-nc-sa/3.0/
% Short Disclaimer: this script is for educational purpose only.
% Longer Disclaimer see  http://igorcarron.googlepages.com/disclaimer
% Modified after suggestions by Ming-Jun Lai, 2017 - 2018
%
% Modified again to use lsqr, as suggested by Tropp and Needell in CoSaMP
% paper. 9 April 2019.

% Initialization
Sest = zeros(size(Phi,2),1);
utrue = Sest;
v = u;
prevresen = norm(v);
t = 1;
numericalprecision = 1e-12;
T2 = [];
while (t <= maxiterations) && (norm(v)/norm(u) > tol)
  y = abs(Phi'*v);
  [vals,~] = sort(y,'descend');
  Omega = find(y >= vals(K) & y > numericalprecision);
  K = length(Omega);
  if K == 0
      Sest = zeros(size(Phi,2),1);
      return
  end 
  T = union(Omega,T2);
  phit = Phi(:,T);
  [bb, ~] = lsqr(phit,u,[],10);  
  b = abs(bb);
  [vals,~] = sort(b,'descend');
  Sest = zeros(length(utrue),1);
  Sest(T(b >= vals(K) & b > numericalprecision)) = b(b >= vals(K) & b > numericalprecision);
  [~,z] = sort(Sest,'descend');
  Told = T2;
  T2 = z(1:K);
  phit = Phi(:,T2);
  [b,~] = lsqr(phit,u,[],10);
  Sest = zeros(length(utrue),1);
  Sest(T2) = b;
  v = u - Phi(:,T2)*b;
  newresen = norm(v);
  if newresen > prevresen
    T2 = Told;
    phit = Phi(:,T2);
    [b,~] = lsqr(phit,u,[],10);
    Sest = zeros(length(utrue),1);
    Sest(T2) = b;
    break;
  end
  prevresen = newresen;
  t = t+1;
end