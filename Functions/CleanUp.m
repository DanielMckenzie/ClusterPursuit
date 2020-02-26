function [ClustersFinal,NewClass,Leftovers] = CleanUp(A,Clusters)
%        ClustersFinal = CleanUp(A,Clusters)
% This function cleans up the clustering found by Cluster Pursuit.
% Attempting to speed it up.
% INPUT
% ===============================
% A .................... Adjacency matrix
% Clusters ............. CELL. Contains clusters found by Cluster Pursuit
%
% OUTPUT
% ==============================
% ClustersFinal ........ Cleaned up clusters
%
% Daniel Mckenzie
% 03 February 2018
%

% ========= Some initialization/parameters ============== %
n = size(A,1);
k = length(Clusters);
AllClassified = Cell2Vec(Clusters);
Leftovers = setdiff(1:n,AllClassified);
m = length(Leftovers);

% ====== Classify Leftovers using a majority-rules scheme ===== % 
Scores = zeros(m,k);
Atemp = A(Leftovers,:);
for a = 1:k
    Scores(:,a) = sum(Atemp(:,Clusters{a}),2);
end

ClustersFinal = cell(size(Clusters));
[~,NewClass] = max(Scores,[],2);

for a = 1:k
    temp = Leftovers(NewClass == a);
    ClustersFinal{a} = union(Clusters{a},temp);
end
end

    