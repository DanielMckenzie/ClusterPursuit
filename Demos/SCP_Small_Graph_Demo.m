% ====================================================================== %
% This script demonstrates SingleClusterPursuit on a small graph
% Daniel Mckenzie 
% 22nd May 2019
% ====================================================================== %
clear, close all, clc
addpath(genpath('../Utilities'),'../Functions')

% ============== Generating adjacency matrix of graph ============== %
n = 150;          % Set the number of vertices of the graph
k = 2;           % number of clusters
n0 = ceil(n/k);   % size of each cluster (equally sized)
p = 0.1;% in-cluster connection probability
q = 0.005;   % between cluster connection probability
A = generateA(n,n0,p,q);
% uncomment below if you wish to visualize the adjacency matrix
%Im1 = mat2gray(full(A));
%imshow(imcomplement(Im1));
%title('The ground truth adjacency matrix')

% ================ Randomly permute the adjacency matrix =============== %
perm = randperm(n);
A = A(perm,perm);
[~,permInv] = sort(perm);
TrueCluster = permInv(1:n0);   % the ground truth first cluster, after permutation.

% ========================== Draw seed vertices =================== %
Gamma = datasample(TrueCluster,5,'Replace',false); % seed vertices

% ========= Visualize the graph with seed vertices highlighted ======== %
G = graph(A);
figure
H = plot(G,'Layout','force','MarkerSize',4);
highlight(H,Gamma,'NodeColor','r','MarkerSize',8);
title('Graph with seed vertices highlighted','FontSize',14)


% ====================== Run SSCP ================================= %
% ====== Parameters
epsilon = 0.2;    % epsilon parameter
reject = 0.4;     % change this to increase/decrease probaility of false negative

tic
Cluster = CP_RWT(A,Gamma,n0,epsilon,3,reject);
time = toc;

% ================= Assess Accuracy ===================== %
accuracy = 100*length(intersect(Cluster,TrueCluster))/n0;
disp(['Found Cluster 1 with an accuracy of ',num2str(accuracy),'%'])

% ===============  Replot graph =========================== %
figure
H2 = plot(G,'Layout','force','MarkerSize',4);
highlight(H2,Cluster,'NodeColor','r','MarkerSize',4);
highlight(H2,Gamma,'NodeColor','r','MarkerSize',8);
title('Graph after Cluster Pursuit. Found Cluster highlighted')
