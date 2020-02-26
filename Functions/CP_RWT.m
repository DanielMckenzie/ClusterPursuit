function Cluster = CP_RWT(A,Gamma,n0,epsilon,t,reject,sigma)
%        Cluster = SingleClusterPursuit(A,Gamma,n0,delta,t,reject)
% SingleClusterPursuit, as described in the Dissertation of Mckenzie.
% Note this implementation uses a random walk to find Omega, and uses
% MATLAB's lsqr in place of backslash in the call to subspacepursuit.
% Daniel Mckenzie
% 13 April 2019

% Renamed to the name convention used in current draft of paper Oct 8th
% 2019

% INPUT
% =================================
% A .................... Adjacency matrix of data converted to graph form.
% Gamma ................ VECTOR. Labelled data within cluster of interest
% n0 ................... (estimated) size of C_a
% epsilon .............. Omega_a will be of size (1+epsilon)n0(a)
% reject ............... value of threshold for ClusterPursuit
% t .................... Depth of random walk 
% sigma ................ sparsity multiplier. Default is 0.13
%
% OUTPUT
% ================================
% Cluster...... VECTOR. The elements in the cluster of interest.
%
if nargin == 6
    sigma = 0.13
end

% ========= Initialization ================= %
n = size(A,1); % number of vertices
degvec = sum(A,2);
Dinv = spdiags(1./degvec,0,n,n);
DinvA = Dinv*A;
L = speye(n,n) - DinvA; 

% ============ Call the subroutines ============= %
Omega = RandomWalkThresh(A,Gamma,n0,epsilon,t);
[Cluster,v] = ClusterPursuit(L,Gamma,Omega,ceil(sigma*n0),reject);
end