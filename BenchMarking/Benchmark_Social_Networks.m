% ==================================================================== %
% Benchmarking Cluster Extraction Algorithms on Social Networks
% Using facebook100 data set and the clusters recommended in:
%   " Capacity Releasing Diffusion for Speed and Locality "
% by Wang, Fountoulakis et al.
% Daniel Mckenzie
% 20th September 2019
% ==================================================================== %

clear, clc, close all
addpath(genpath('../ThirdParty'),genpath('../Utilities'),'../Functions')
addpath('../../facebook100')

% ============================= Parameters ===================== %
num_trials = 50;

for schools = 3
    if schools == 1
        load('JohnsHopkins_with_TrueClusters.mat')
    elseif schools == 2
        load('Rice_with_TrueClusters.mat')
    elseif schools == 3
        load('Simmons_with_TrueClusters.mat')
    elseif schools == 4
        load('Colgate_with_TrueClusters.mat')
    end
    n = size(A,1);
    
    % ============== Define all vectors of interest =========== %
    time_RWT_vec = zeros(num_trials,1);
    time_CP_vec = zeros(num_trials,1);
    time_HKGrow_vec = zeros(num_trials,1);
    time_LBSA_vec = zeros(num_trials,1);
    time_ESSC_vec = zeros(num_trials,1);
    time_PPR_vec = zeros(num_trials,1);

    Jaccard_RWT_vec = zeros(num_trials,1);
    Jaccard_CP_vec = zeros(num_trials,1);
    Jaccard_HKGrow_vec = zeros(num_trials,1);
    Jaccard_ESSC_vec = zeros(num_trials,1);
    Jaccard_LBSA_vec = zeros(num_trials,1);
    Jaccard_PPR_vec = zeros(num_trials,1);
    
    Precision_RWT_vec = zeros(num_trials,1);
    Precision_CP_vec = zeros(num_trials,1);
    Precision_HKGrow_vec = zeros(num_trials,1);
    Precision_ESSC_vec = zeros(num_trials,1);
    Precision_LBSA_vec = zeros(num_trials,1);
    Precision_PPR_vec = zeros(num_trials,1);
    
    Recall_RWT_vec = zeros(num_trials,1);
    Recall_CP_vec = zeros(num_trials,1);
    Recall_HKGrow_vec = zeros(num_trials,1);
    Recall_ESSC_vec = zeros(num_trials,1);
    Recall_LBSA_vec = zeros(num_trials,1);
    Recall_PPR_vec = zeros(num_trials,1);
    
    if schools == 1
        TrueCluster = TrueClusters{2};
    else 
        TrueCluster = TrueClusters{2};
    end
    n0 = length(TrueCluster);
    
    for j = 1:num_trials
         % ================ Draw Seed sets =============== %
         Gamma = datasample(TrueCluster,ceil(0.02*n0),'Replace',false);
         Degs_in_cluster = sum(A(TrueCluster,:),2);
         [~,Gamma1] = max(Degs_in_cluster); 
         Gam_NBD = find(A(Gamma1,:));
         
         % ============ Find Cluster with only RWThresh =========== %
         tic
         Cluster_RWT = RandomWalkThresh(A,Gamma,n0,0,3);
         time_RWT_vec(j) = toc;
         Jaccard_RWT_vec(j) = Jaccard_Score(TrueCluster,Cluster_RWT);
         Precision_RWT_vec(j) = Precision(TrueCluster,Cluster_RWT);
         Recall_RWT_vec(j) = Recall(TrueCluster,Cluster_RWT);
         
         % ========== Find Cluster with ClusterPursuit ============ %
        tic
        Cluster_CP = CP_RWT(A,Gamma,n0,0.25,3,0.5,0.5);
        time_CP_vec(j) = toc;
        Jaccard_CP_vec(j) = Jaccard_Score(TrueCluster,Cluster_CP);
        Precision_CP_vec(j) = Precision(TrueCluster,Cluster_CP);
        Recall_CP_vec(j) = Recall(TrueCluster,Cluster_CP);

        % ============= ESSC ==================== %
%          tic
%          Cluster_ESSC  = ESSC(A,Gam_NBD,0.01);
%          time_ESSC_vec(j) = toc
%          Jaccard_ESSC_vec(j) = Jaccard_Score(TrueCluster,Cluster_ESSC)
%          Precision_ESSC_vec(j) = Precision(TrueCluster,Cluster_ESSC);
%          Recall_ESSC_vec(j) = Recall(TrueCluster,Cluster_ESSC);

        % ================== HKGrow =================== %
        tic
        [Cluster_HKGrow,~,~,~] = hkgrow(A,Gamma);
        time_HKGrow_vec(j) =  toc;
        Jaccard_HKGrow_vec(j) = Jaccard_Score(TrueCluster,Cluster_HKGrow);
        Precision_HKGrow_vec(j) = Precision(TrueCluster,Cluster_HKGrow);
        Recall_HKGrow_vec(j) = Recall(TrueCluster,Cluster_HKGrow);

        % ========== Find Cluster with LBSA algorithm ========== %
        tic
        Cluster_LBSA = LBSA2(A,Gamma,'hk','Lanczos');
        time_LBSA_vec(j) = toc;
        Jaccard_LBSA_vec(j) = Jaccard_Score(TrueCluster,Cluster_LBSA);
        Precision_LBSA_vec(j) = Precision(TrueCluster,Cluster_LBSA);
        Recall_LBSA_vec(j) = Recall(TrueCluster,Cluster_LBSA);
        
        % ============= Personal Page Rank ===================== %
        tic
        expectedVolMult = 90*n0/length(Gamma);  % Volume multiple. Calculated as av_degree*size_cluster/num_seed_vertices
        degs = sum(A,1);
        Dinv = diag(degs.^(-1));
        L = eye(n,n) - Dinv*A;
        [V,D] = eigs(L,2,'smallestreal');
        alpha_param = D(2,2);
        [Cluster_PPR,~,~,~] = pprgrow(A,Gamma,'alpha',4*alpha_param);
        time_PPR_vec(j) = toc;
        Jaccard_PPR_vec(j) = Jaccard_Score(TrueCluster,Cluster_PPR);
        Precision_PPR_vec(j) = Precision(TrueCluster,Cluster_PPR);
        Recall_PPR_vec(j) = Recall(TrueCluster,Cluster_PPR);
    end
    
%     figure
%     Time = [time_RWT_vec,time_CP_vec,time_HKGrow_vec,time_LBSA_vec,time_PPR_vec,time_ESSC_vec];
%     boxplot(log(Time),'Labels',{'RWT','CP+RWT','HKGrow','LBSA','PPR','ESSC'})
%     ylabel('Log. of Run Time','FontSize',14)
%     set(gca,'FontSize',14)
%     if schools == 1
%         title('Johns Hopkins: Class of 2009')
%     elseif schools == 2
%         title('Rice: Dorm 203')
%     elseif schools == 3
%         title('Simmons College: Class of 2009')
%     else
%         title('Colgate University: Class of 2006')
%     end
%     
%     figure
%     Jaccard = [Jaccard_RWT_vec,Jaccard_CP_vec,Jaccard_HKGrow_vec,Jaccard_LBSA_vec,Jaccard_PPR_vec,Jaccard_ESSC_vec];
%     boxplot(Jaccard,'Labels',{'RWT','CP+RWT','HKGrow','LBSA','PPR','ESSC'})
%     %ylabel('Jaccard Index','FontSize',14)
%     set(gca,'FontSize',14)
%     if schools == 1
%         title('Johns Hopkins: Class of 2009')
%     elseif schools == 2
%         title('Rice: Dorm 203')
%     elseif schools == 3
%         title('Simmons College: Class of 2009')
%     else
%         title('Colgate University: Class of 2006')
%     end
end

% figure
% Precision = [Precision_RWT_vec,Precision_CP_vec,Precision_HKGrow_vec,Precision_LBSA_vec,Precision_PPR_vec,Precision_ESSC_vec];
% boxplot(Precision,'Labels',{'RWT','CP+RWT','HKGrow','LBSA','PPR','ESSC'})
% ylabel('Precision','FontSize',14)
% 
% 
% figure
% Recall = [Recall_RWT_vec,Recall_CP_vec,Recall_HKGrow_vec,Recall_LBSA_vec,Recall_PPR_vec,Recall_ESSC_vec];
% boxplot(Recall,'Labels',{'RWT','CP+RWT','HKGrow','LBSA','PPR','ESSC'})
% ylabel('Recall','FontSize',14)

figure 
Precision_and_Recall = [Precision_CP_vec,Recall_CP_vec,Precision_HKGrow_vec,Recall_HKGrow_vec,Precision_LBSA_vec,Recall_LBSA_vec,Precision_PPR_vec,Recall_PPR_vec];
boxplot(Precision_and_Recall,'Labels',{'CP+RWT:Prec.','CP+RWT:Rec.','HKGrow:Prec.','HKGrow:Rec.','LBSA:Prec.','LBSA:Rec.','PPR:Prec.','PPR:Rec.'})
set(gca,'FontSize',10,'XTickLabelRotation',90)