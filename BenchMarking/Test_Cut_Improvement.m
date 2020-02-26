% ==================================================================== %
% Test Cut Improvement
% Daniel Mckenzie
% 4th October 2019
% ==================================================================== %

clear, clc, close all
addpath(genpath('../ThirdParty'),genpath('../Utilities'),'../Functions')

% ==================== Parameters ==================== %
% ===== for Model
n1 = 500;
n0vec = n1*[1,1.5,2,2.5];
n = sum(n0vec);
p_prime = ((log(n))^2)/5;
q = log(n)/n;
P_diag = p_prime./n0vec - q;
P = q*ones(4,4) + diag(P_diag);

% ==== For RWThresh
sample_frac = 0.02;
epsilon = 0.13;

% ==== For ClusterPursuit
reject = 0.5;

% ==== For SimpleLocal


% ===================== Generate Graph ==================== %
A = generateA2(n0vec,P);
%Im1 = mat2gray(full(A));
perm = 1:n;
perm = randperm(n);
A = A(perm,perm);
[~,permInv] = sort(perm);
TrueCluster = permInv(1:n1);
Gamma = datasample(TrueCluster,ceil(sample_frac*n1),'Replace',false);


% ==================== Run Algorithms ================ %
% ==== Generate initial cut with RWThresh.
Omega = RandomWalkThresh(A,Gamma,n1,2*epsilon,3);
RWT_Jac = Jaccard_Score(TrueCluster,Omega)

% ==== Use ClusterPursuit
degvec = sum(A,2);
Dinv = spdiags(1./degvec,0,n,n);
DinvA = Dinv*A;
L = speye(n,n) - DinvA; 
[C1_SCP,v] = ClusterPursuit2(L,Gamma,Omega,n1,reject);

SCP_Jac = Jaccard_Score(TrueCluster,C1_SCP)

% === Use SimpleLocal
[C1_SimpleLocal_vec,~] = SimpleLocal(A,Omega,0.5);
C1_SimpleLocal = find(C1_SimpleLocal_vec);

SimpleLocal_Jac = Jaccard_Score(TrueCluster,C1_SimpleLocal)


% === Use Flow Improve
Omega_indicator = false(n,1);
Omega_indicator(Omega) = true;

[C1_FlowImprove_vec,~] = flow_improve(A,Omega_indicator);
C1_FlowImprove = find(C1_FlowImprove_vec);

FlowImprove_Jac = Jaccard_Score(TrueCluster,C1_FlowImprove)



