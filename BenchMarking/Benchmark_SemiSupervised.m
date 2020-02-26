% ===================================================================== %
% Testing Semi-Supervised classification with SingleClusterPursuit on MNIST
% Can use KNN_Max adjacecny matrix as well as shortest paths based
% adjacency matrices.
% Daniel Mckenzie
% 8 April 2019
% ===================================================================== %

clear, clc, close all
addpath('../Functions',genpath('../Utilities'))
addpath('../Mat_files')

% =================== Parameters and load the data ===================== %
%load('MNIST_KNN_Max_K=15.mat');
%load('MNIST_KNN_Mult_K=15.mat');
load('OptDigits_KNN_Mult_K=15.mat');
%load('COIL_Matrices')
%A = AMult;
%load('ShortestPaths_Adj_Matrices2.mat')
%A = A1;
%clear A1 A2 A3 A4
%y = labels;
%y = ynew;

% === for ICP+RWT
epsilon = 0.13;
reject = 0.5;
num_trials = 10;
num_sizes = 5;

% === For three moons
%r1 = 1;
%r2 = 1;
%r3 = 1.5;
%ambient_dim = 100;
%noise_level = 0.14;
%n0 = 500;

k = 10;  %number of clusters
n = size(A,1);
%n = 1500;

% =========== Find the ground truth clusters ======== %
TrueClusters = cell(k,1);
%n0vec = n0*ones(3,1);
 for a = 1:k
     Ctemp = find(y== a-1);
     TrueClusters{a} = Ctemp;
     n0vec(a) = length(Ctemp);
     
 end

% ======== Define all Vectors of Interest ============= %
time_ICP_vec = zeros(num_sizes,1);
Acc_ICP_vec = zeros(num_sizes,1);

for j = 1:num_sizes
    sample_frac = 0.005*(j);
    time_ICP = 0;
    Acc_ICP = 0;
    for i = 1:num_trials
        %perm = randperm(n);
        %[~,permInv] = sort(perm);
        %TrueClusters{1} = permInv(1:500);
        %TrueClusters{2} = permInv(501:1000);
        %TrueClusters{3} = permInv(1001:1500);
        %[~,Points] = Generate3Moons(r1,r2,r3,n0,noise_level,ambient_dim);
        %A = CreateKNN_Mult_from_Data(Points,15,10);
        %A = A(perm,perm);
        Gamma = cell(k,1);
        % Randomly sample labeled data
        for a = 1:k
            Gamma{a} = datasample(TrueClusters{a},ceil(sample_frac*n0vec(a)),'Replace',false);
        end
        tic
        [~,Clusters_ICP] = IteratedClusterPursuit(A,Gamma,n0vec,epsilon,reject,5);
        time_ICP = time_ICP + toc;
        
        Corr_Class_ICP = 0;
        for a = 1:k
            Corr_Class_ICP = Corr_Class_ICP + length(intersect(Clusters_ICP{a},TrueClusters{a}));
        end
        Acc_ICP = Acc_ICP + (Corr_Class_ICP/n)
    end
    time_ICP_vec(j) = time_ICP/num_trials;
    Acc_ICP_vec(j) = Acc_ICP/num_trials;
end

    