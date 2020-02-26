function [Clusters,ClustersFinal] = IteratedClusterPursuit(A,Gamma,n0vec,epsilon,reject,t)
%        Clusters = ISSCP2(A,Gamma,n0vec,delta,reject,weight))
% This function iterates the SSCP function to find clusters of sizes
% specified by n0vec. It does not allow for overlaps as it removes the 
% i-th cluster before attempting to find the (i+1)-st cluster. 
%
% This version IS PARTITIONAL. That is, every vertex is assigned to
% precisely one cluster.
% Experimental version---uses ClusterPursuit2
% Daniel Mckenzie 
% September 20th 2019
%
% INPUT
% =================================
% A .................... Adjacency matrix of data converted to graph form.
% Gamma ......... CELL. Gamma{a} is labelled data within C_a
% n0vec ............ VECTOR. n0vec(a) is the (estimated) size of C_a
% epsilon ......... Omega_a will be of size (1+epsilon)n0vec(a)
% reject ....... value of threshold for ClusterPursuit
% t ........... Depth of Random Walk
%
% OUTPUT
% ================================
% Clusters ...... CELL. Clusters{a} = C_a before applying CleanUp
% ClustersFinal ....... CELL. ClustersFinal{a} = C_a after applying
% CleanUp.

% Note that the CleanUp sub-routine corresponds is not explicitly discussed
% in the paper. Essentially it just uses a k-nearest neighbor classifier to
% assign the ``leftover'' vertices, i.e. those in V\C_1U...UC_k. If this is
% not desirable one can use Clusters instead of ClustersFinal.
%
% Daniel Mckenzie
% 02 February 2018
% Updated 14 April 2019 to use lsqr instead of backslash in call to
% subspacepursuit.
%

% ========= Initialization ================= %
Aold = A;
n = size(A,1); % number of vertices
N = n; % number of vertices, immutable
k = length(n0vec); %number of clusters
Clusters = cell(k,1);
Remaining = cell(k,1);
OmegaCell = cell(k,1);

sigma = 0.26;  % Hard code the sparsity parameter. Can change this if necessary.

% ============ Now iteratively find the clusters ============ %
for a=1:k-1
    a
    % ====== Now run through the main SCP algorithm ====== %
    C = CP_RWT(A,Gamma{a},n0vec(a),epsilon,t,reject,sigma);
    Remaining{a} = setdiff(1:n,C);
    
    if a ~=1
        for b = a-1:-1:1 
            Remtemp = Remaining{b};
            Ctemp = Remtemp(C);
            C = Ctemp;
        end
    end
    Clusters{a} = C;
    % also need to readjust the indices of the labeled sets Gamma_b
    for b = a+1:k-1
        Gamtemp = find(ismember(Remaining{a},Gamma{b}));
        Gamma{b} = Gamtemp';
    end
    if a <= k-1
        A = A(Remaining{a},Remaining{a});
        n = size(A,1);
    end
end

% ====================== Extract the final cluster ===================== %
Dinv = spdiags(1./sum(Aold,2),0,N,N);
DinvA = Dinv*Aold;
L = speye(N,N) - DinvA;
Omega = setdiff(1:N,Cell2Vec(Clusters(1:k-1)));
[C,~] = ClusterPursuit(L,Gamma{k},Omega,n0vec(k),reject);
C = setdiff(C,Cell2Vec(Clusters(1:k-1)));  % remove any previously classified vertices that may have snuck back in 
Clusters{k} = C;

% ================ Do a final sweep and classify the leftovers ========== %
[ClustersFinal,~,~] = CleanUp(Aold,Clusters);

end

    
            
            
        
    
