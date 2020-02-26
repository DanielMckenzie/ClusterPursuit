function [C,v] = ClusterPursuit(L,Gamma_a,Omega_a,s,reject)
%        [C, v] = ClusterPursuit(L,Gamma_a,Omega_a,n_a)
% Subroutine for SingleClusterPursuit that performs the pursuit step.
% Uses lsqr in place of backslash in the call to subspacepursuit.
% 
% Experimental version 19th September --- allowing algorithm to add
% vertices that are not in Omega.
%
% INPUT
% ==========================================
% L .................... Laplacian matrix
% Gamma_a .............. Labeled data for C_a
% Omega_a .............. Superset containing C_a
% s .................... Sparsity - estimated size of C_a/Omega + Omega/C_a
% reject ............... number in (0,1), value at which we threshold.
%
% OUTPUT
% =========================================
% C ................... Estimate of C_a (including Gamma_a)
% v ................... Output of subspacepursuit - vector of probabilities
% of not being in C
%

%Phi = L(:,Omega_a);
Phi = L;
n = size(L,1);
yg = sum(Phi(:,Omega_a),2);
g = length(Gamma_a);
%sparsity = ceil(2*(length(Omega_a) - n_a - g));
sparsity = s;
if sparsity <= 0
    C = union(Omega_a,Gamma_a);
    v = zeros(n,1);
else
    v = subspacepursuit(Phi,yg,sparsity,1e-10,ceil(log(10*n)));
    W_a = find(v > reject);
    U_a  = find(v <- reject);
    C = setdiff(Omega_a,W_a);
    C = union(C,Gamma_a);
    C = union(C,U_a);
end

end
