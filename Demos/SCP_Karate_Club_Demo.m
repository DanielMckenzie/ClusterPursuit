% ===================================================================== %
% This script demonstrates the SingleClusterPursuit algorithm on the 
% commonly used Zachary's karate club data set.
% Daniel Mckenzie
% 22nd May 2019
% ===================================================================== %

clear, close all, clc
addpath(genpath('../Utilities'),'../Functions','../Mat_files')
load('Karate.mat')

% ======================== Choose the seed vertices ==================== %
% Different choices of Gamma can lead to very different results. Below are
% a few interesting options. Note that one can also vary the parameters
% epsilon, reject and n0 to get different results. 
Gamma = [30,33,24];
%Gamma = [30,33,34]
% Gamma = [1,8,11];

% ========= Visualize the graph with seed vertices highlighted ======== %
G = graph(A);
figure
H = plot(G,'Layout','force','MarkerSize',4);
return
highlight(H,Gamma,'NodeColor','r','MarkerSize',8);
title('Zacharys Karate Club, with three seed vertices in red.','FontSize',14)

% ======================= Ground Truth ============================ %
TrueCluster1 =[1,2,3,4,5,6,7,8,11,12,13,14,17,18,20,22];
TrueCluster2 = setdiff(1:34,TrueCluster1);

% ====================== Run SSCP ================================= %
% ====== Parameters
epsilon = 0.2;    % epsilon parameter
reject = 0.8;     % change this to increase/decrease probaility of false negative
n0 = 17;          % target cluster size


tic
Cluster = CP_RWT(A,Gamma,n0,epsilon,2,reject);
time = toc;

% ====================== Determine accuracy =================== %
Jaccard = Jaccard_Score(Cluster,TrueCluster2);
disp(['Found Cluster 1 with an accuracy of ',num2str(100*Jaccard),'%'])

% ===============  Replot graph =========================== %
figure
H2 = plot(G,'Layout','force','MarkerSize',4);
highlight(H2,Cluster,'NodeColor','r','MarkerSize',4);
highlight(H2,Gamma,'NodeColor','r','MarkerSize',8);
title('Zacharys Karate Club, with cluster containing seed vertices in red.')
