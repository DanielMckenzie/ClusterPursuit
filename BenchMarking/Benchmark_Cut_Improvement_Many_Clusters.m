% ==================================================================== %
% Benchmarking two venerable cut improvement algorithms: FlowImprove and
% SimpleLocal against ClusterPursuit.
% Daniel Mckenzie
% October 8th 2019
% ==================================================================== %

clear, clc, close all
addpath(genpath('../ThirdParty'),genpath('../Utilities'),'../Functions')

% ==================== Parameters ==================== %
num_trials = 10;
num_sizes = 5;
Cluster_sizes = 200+100*[1:num_sizes];  % Vector of cluster sizes

% ==== For RWThresh
sample_frac = 0.01;
epsilon = 0.13;

% ==== For ClusterPursuit
reject = 0.5;

% ==== For SimpleLocal

% ============ Create matrices of interest ================= %
time_CP_mat = zeros(num_trials,num_sizes);
time_SimpleLocal_mat = zeros(num_trials,num_sizes);
time_FlowImprove_mat = zeros(num_trials,num_sizes);

Jaccard_RWT_mat = zeros(num_trials,num_sizes);
Jaccard_CP_mat = zeros(num_trials,num_sizes);
Jaccard_SimpleLocal_mat = zeros(num_trials,num_sizes);
Jaccard_FlowImprove_mat = zeros(num_trials,num_sizes);



for j = 1:num_sizes
    n1 = Cluster_sizes(j);
    n0vec = n1*[1,1.5,2,2.5];
    n = sum(n0vec);
    p_prime = ((log(n))^2)/2;
    q = 5*log(n)/n;
    P_diag = p_prime./n0vec - q;
    P = q*ones(4,4) + diag(P_diag);
    
    %n0vec = [n1,9*n1];
    %n = sum(n0vec);
    %q = log(n)/n;
    %P = [(log(n1)^2)/n1,q;q,q];
    

    for i = 1:num_trials
        % ========= Generate graph
        A = generateA2(n0vec,P);
        perm = randperm(n);
        A = A(perm,perm);
        [~,permInv] = sort(perm);
        TrueCluster = permInv(1:n1);
        Gamma = datasample(TrueCluster,ceil(sample_frac*n1),'Replace',false);

        % ==================== Run Algorithms ================ %
        % ==== Generate initial cut with RWThresh.
        Omega = RandomWalkThresh(A,Gamma,n1,epsilon,3);
        Jaccard_RWT_mat(i,j) = Jaccard_Score(TrueCluster,Omega);

        % ==== Use ClusterPursuit
        degvec = sum(A,2);
        Dinv = spdiags(1./degvec,0,n,n);
        DinvA = Dinv*A;
        L = speye(n,n) - DinvA;
        tic
        [C1_CP,~] = ClusterPursuit(L,Gamma,Omega,ceil(0.13*n1),reject);
        time_CP_mat(i,j) = toc;

        Jaccard_CP_mat(i,j) = Jaccard_Score(TrueCluster,C1_CP);

        % === Use SimpleLocal
        tic
        [C1_SimpleLocal_vec,~] = SimpleLocal(A,Omega,0.5);
        time_SimpleLocal_mat(i,j) = toc;
        
        C1_SimpleLocal = find(C1_SimpleLocal_vec);

        Jaccard_SimpleLocal_mat(i,j) = Jaccard_Score(TrueCluster,C1_SimpleLocal);


        % === Use Flow Improve
        Omega_indicator = false(n,1);
        Omega_indicator(Omega) = true;

        tic
        [C1_FlowImprove_vec,~] = flow_improve(A,Omega_indicator);
        time_FlowImprove_mat(i,j) = toc;
        
        C1_FlowImprove = find(C1_FlowImprove_vec);
        Jaccard_FlowImprove_mat(i,j) = Jaccard_Score(TrueCluster,C1_FlowImprove);
        
    end
end

% ======= Plot all for comparison ======== %
x_vals = Cluster_sizes;
figure, hold on
%plot(x_vals,mean(Jaccard_SimpleLocal_mat,1),'LineWidth',3)
plot(x_vals,mean(Jaccard_FlowImprove_mat,1),'LineWidth',3)
plot(x_vals,mean(Jaccard_RWT_mat,1),'LineWidth',3)
plot(x_vals,mean(Jaccard_CP_mat,1),'k','LineWidth',3)

legend({'FlowImprove','RWT','ClusterPursuit'},'FontSize',14)
ylabel('Jaccard Index')
xlabel('Size of n_1')
set(gca, 'FontSize',14)

% ======= Plot all for times comparison ======== %
figure, hold on
%plot(x_vals,log(mean(time_SimpleLocal_mat,1)),'LineWidth',3)
plot(x_vals,log(mean(time_FlowImprove_mat,1)),'LineWidth',3)
plot(x_vals,log(mean(time_CP_mat,1)),'k','LineWidth',3)

legend({'FlowImprove','ClusterPursuit'},'FontSize',14)
ylabel('logarithm of run time')
xlabel('Size of n_1')
set(gca, 'FontSize',14)

