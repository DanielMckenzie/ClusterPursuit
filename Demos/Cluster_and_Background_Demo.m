% ========================= %
% This script demonstrates  how to use SingleClusterPursuit to extract a 
% single cluster from a graph drawn at random from the Stochastic Block
% Model with one cluster and a large background.
% Daniel Mckenzie and Ming-Jun Lai
% 22nd May 2019
% ======================== %

%clear, close all, clc
addpath(genpath('../Utilities'),'../Functions')

% =============== Parameters ================= %
n0 = 1000;
n0vec = [n0,5000];
k = 2;
epsilon = 0.2;
reject = 0.8;
n = sum(n0vec);
p = 5*(log(n))^2/(n);% in-cluster connection probability
q = log(n)/n;   % between cluster connection probability
P = [p,q;q,q];

% ============= Now generate matrix and display in greyscale ============ %
A = generateA2(n0vec,P);
Im1 = mat2gray(full(A));
imshow(imcomplement(Im1));
title('The ground truth adjacency matrix')

% ================ Randomly permute the adjacency matrix =============== %
perm = randperm(n);
A = A(perm,perm);
[~,permInv] = sort(perm);
TrueCluster = permInv(1:n0);   % the ground truth first cluster, after permutation.
Im2 = mat2gray(full(A));
figure
imshow(imcomplement(Im2));
title('The adjacency matrix, randomly permuted')

% =============== Run SingleClusterPursuit  ======================== %
Gamma = datasample(TrueCluster,5,'Replace',false); % the set of seed vertices

tic
Cluster = CP_RWT(A,Gamma,n0,epsilon,3,reject);
time = toc;

accuracy = 100*length(intersect(Cluster,TrueCluster))/n0;

NewInds = [Cluster', setdiff(1:n,Cluster)];
Anew = A(NewInds,NewInds);
Im3 = mat2gray(full(Anew));
figure
imshow(imcomplement(Im3))
title('The adjacency matrix, permuted to reveal the cluster found')

disp(['Found Cluster 1 in ', num2str(time), ' seconds, with ', num2str(accuracy), '% accuracy'])
