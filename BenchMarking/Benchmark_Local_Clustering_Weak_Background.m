% ================================================================= %
% Test Experiment on SBM. One cluster and weak background
% Daniel Mckenzie
% October 2019
% ================================================================= %

clear, clc, close all
addpath(genpath('../ThirdParty'),genpath('../Utilities'),'../Functions')

% ============== Parameters ================= %
num_sizes = 5;                      % Number of different cluster sizes
num_trials = 5;                    % Number of trials to run for each size
Cluster_sizes = 200+100*[1:num_sizes];  % Vector of cluster sizes


% ============ Parameters for various algorithms ========== %
% === for CP+RWT
epsilon = 0.1;
reject = 0.5;
sample_frac = 0.01;  % fraction of True Cluster to sample 
% === for LOSP_Plus
alpha = 0.1;         % random walk diffusion parameter
WalkMode = 2;        % type of random walk for LOSP++
d = 2;               % dimension of local spectral subspace
kk = 2;              % number of random walk steps

% ============== Define all matrices of interest =========== %
time_RWT_mat = zeros(num_trials,num_sizes);
time_CP_mat = zeros(num_trials,num_sizes);
time_HKGrow_mat = zeros(num_trials,num_sizes);
time_LBSA_mat = zeros(num_trials,num_sizes);
time_ESSC_mat = zeros(num_trials,num_sizes);
time_PPR_mat = zeros(num_trials,num_sizes);

Jaccard_RWT_mat = zeros(num_trials,num_sizes);
Jaccard_CP_mat = zeros(num_trials,num_sizes);
Jaccard_HKGrow_mat = zeros(num_trials,num_sizes);
Jaccard_LBSA_mat = zeros(num_trials,num_sizes);
Jaccard_ESSC_mat = zeros(num_trials,num_sizes);
Jaccard_PPR_mat = zeros(num_trials,num_sizes);

for j = 1:num_sizes
    n1 = Cluster_sizes(j);
    n0vec = n1*[1,10];
    n = sum(n0vec);
    p_prime = ((log(n))^2)/2;
    q = 10*log(n);
    P = [p_prime/n1, q/n; q/n, q/n];
    
    for i = 1:num_trials
        A = generateA2(n0vec,P);
        %Im1 = mat2gray(full(A));
        perm = randperm(n);
        A = A(perm,perm);
        
        % =============== Find ground truth Cluster ================ %
        [~,permInv] = sort(perm);
        TrueCluster = permInv(1:n1);
        
        % ============== ExtractSeed vertices ================ %
        Gamma1 = datasample(TrueCluster,1,'Replace',false);
        Gamma = datasample(TrueCluster,ceil(sample_frac*n1),'Replace',false);
        Gam_NBD = find(A(Gamma1,:));
        
        % ========== Find Cluster with only RandomWalkThresh =========== %
        tic 
        Cluster_RWT = RandomWalkThresh(A,Gamma,n1,0,3);
        time_RWT_mat(i,j) = toc;
        Jaccard_RWT_mat(i,j) = Jaccard_Score(TrueCluster,Cluster_RWT)
        
        % ========== Find Cluster with ClusterPursuit ============ %
        tic
        Cluster_CP = CP_RWT(A,Gamma,n1,epsilon,3,reject);
        time_CP_mat(i,j) = toc;
        Jaccard_CP_mat(i,j) = Jaccard_Score(TrueCluster,Cluster_CP)

              
        % ========== Find Cluster with ESSC algorithm ============ %
%          if 1 %n0 <= 200
%              tic
%              Cluster_ESSC  = ESSC(A,Gam_NBD,0.05);
%              time_ESSC_mat(i,j) = toc
%              Jaccard_ESSC_mat(i,j) = Jaccard_Score(TrueCluster,Cluster_ESSC)
%          end
        
        % ========== Find Cluster with HKGrow algorithm ========= %
        tic
        [Cluster_HKGrow,~,~,~] = hkgrow(A,Gamma);
        time_HKGrow_mat(i,j) = toc;
        Jaccard_HKGrow_mat(i,j) = Jaccard_Score(TrueCluster,Cluster_HKGrow)
        
        % ========== Find Cluster with Personalized Page Rank algorithm ========= %
        tic
        degs = sum(A,1);
        Dinv = diag(degs.^(-1));
        L = eye(n,n) - Dinv*A;
        [V,D] = eigs(L,2,'smallestreal');
        alpha_param = D(2,2);
        %expectedVolMult = n1*p_prime/length(Gamma);
        [Cluster_PPR,~,~,~] = pprgrow(A,Gamma,'alpha',alpha_param);
        time_PPR_mat(i,j) = toc;
        Jaccard_PPR_mat(i,j) = Jaccard_Score(TrueCluster,Cluster_PPR)
        
         % ========== Find Cluster with LBSA algorithm ========== %
         if j <= 4
            tic
            Cluster_LBSA = LBSA2(A,Gamma,'rw','Lanczos');
            time_LBSA_mat(i,j) = toc;
            Jaccard_LBSA_mat(i,j) = Jaccard_Score(TrueCluster,Cluster_LBSA)
         end
    end
end

boxplot(Jaccard_CP_mat,'Labels',{'300','400','500','600','700'})
set(gca, 'FontSize',14)

% ======= Plot all for comparison ======== %
figure, hold on
plot([300:100:700],mean(Jaccard_HKGrow_mat,1),'LineWidth',3)
plot([300:100:700],mean(Jaccard_PPR_mat,1),'LineWidth',3)
plot([300:100:700],mean(Jaccard_RWT_mat,1),'LineWidth',3)
plot([300:100:600],mean(Jaccard_LBSA_mat(:,1:4),1),'LineWidth',3)
%plot([300:100:700],mean(Jaccard_ESSC_mat,1),'LineWidth',3)
plot([300:100:700],mean(Jaccard_CP_mat,1),'k','LineWidth',3)
legend({'HKGrow','PPR','RWT','LBSA','CP+RWT'},'FontSize',14)
ylabel('Jaccard Index')
xlabel('Size of n_1')
set(gca, 'FontSize',14)

% ======= Plot all for times comparison ======== %
figure, hold on
plot([300:100:700],log(mean(time_HKGrow_mat,1)),'LineWidth',3)
plot([300:100:700],log(mean(time_PPR_mat,1)),'LineWidth',3)
plot([300:100:700],log(mean(time_RWT_mat,1)),'LineWidth',3)
plot([300:100:600],log(mean(time_LBSA_mat(:,1:4),1)),'LineWidth',3)
%plot([300:100:700],log(mean(time_ESSC_mat,1)),'LineWidth',3)
plot([300:100:700],log(mean(time_CP_mat,1)),'k','LineWidth',3)
legend({'HKGrow','PPR','RWT','LBSA','CP+RWT'},'FontSize',14)
ylabel('logarithm of run time')
xlabel('Size of n_1')
set(gca, 'FontSize',14)

        
        
    