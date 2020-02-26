function Omega = RandomWalkThresh(A,Gamma,n0_hat,epsilon,t)
%        Omega = RandomWalkThresh(A,Gamma,n0_hat,epsilon,t)
% This function determines the superset Omega needed for SSCP by running a
% random walk for k steps starting from Gamma and returning the
% (1+epsilon)n0_hat vertices with highest probabilities of being visited.
% 
% INPUTS
% ==================================================================== %
% A ............ Adjacency matrix
% Gamma ........ Seed vertices/ labelled data
% n0_hat ....... Estimate on the size of the cluster
% epsilon ...... Oversampling parameter
% t ............ Depth of random walk
% 
% OUTPUTS
% ==================================================================== %
% Omega ....... Superset containing a large fraction of the vertices in
% cluster.
% 
% Daniel Mckenzie
% 2 March 2019

% =========================== Initialization ========================== %
n = size(A,1);
Dtemp = sum(A,2);
Dinv = spdiags(Dtemp.^(-1),0,n,n);
v0 = sparse(Gamma,1,Dtemp(Gamma),n,1);
P = A*Dinv;

% ===================== Random Walk and Threshold ===================== %
v = v0;
for i = 1:t
    v = P*v;
end
size(v)
[w,IndsThresh] = sort(v,'descend');
FirstZero = find(w==0, 1, 'first');
if ~isempty(FirstZero) && FirstZero < ceil((1+epsilon)*(n0_hat))
    warning('the size of Omega is smaller than (1+delta) times the user specified cluster size. Try a larger value of k')
    T = FirstZero;
else
    T = ceil((1+epsilon)*(n0_hat));
end
Omega = union(IndsThresh(1:T),Gamma);

