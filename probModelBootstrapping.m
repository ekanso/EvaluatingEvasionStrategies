clear
close all
clc

%% Load Data

load('evasionData.mat');
%% Sort Data
phi = data.Phi(isfinite(data.U));
u   = data.U(isfinite(data.U));
v   = data.V(isfinite(data.U))*0.01;
lam = data.Lambda(isfinite(data.U));
theta = data.Delta(isfinite(data.U));

phi_02 = phi(v == 0.02);
phi_11 = phi(v == 0.11);
phi_20 = phi(v == 0.2);
lam_02 = lam(v == 0.02);
lam_11 = lam(v == 0.11);
lam_20 = lam(v == 0.2);
u_02 = u(v == 0.02);
u_11 = u(v == 0.11);
u_20 = u(v == 0.2);
v_02 = v(v == 0.02);
v_11 = v(v == 0.11);
v_20 = v(v == 0.2);
theta_02 = theta(v == 0.02);
theta_11 = theta(v == 0.11);
theta_20 = theta(v == 0.2);

%% Setup Data

Ns_02 = length(phi_02);

x_02 = [phi_02;...
    lam_02;...
    u_02;...
    v_02;...
    theta_02];
Ns_11 = length(phi_11);

x_11 = [phi_11;...
    lam_11;...
    u_11;...
    v_11;...
    theta_11];
Ns_20 = length(phi_20);

x_20 = [phi_20;...
    lam_20;...
    u_20;...
    v_20;...
    theta_20];

Ns = length(phi);

x = [phi;
    lam;
    u;
    v;
    theta];
%% Run Optimization

options = optimoptions('fmincon','Display','notify-detailed','MaxFunctionEvaluations',1e4);
%% no bootstrapping
% eta0_CL = eta0_s1;
% eta0_AP = eta0_s2;
% eta0_PL = eta0_s4;
% eta0_DO = eta0_s4;
% eta0_OT = eta0_s3;
% eta0_CL = [1;1];
% eta0_AP = [1;1];
% eta0_PL = [1];
% eta0_DO = [1;1];
% eta0_OT = [1;1];
% eta0_CL_02 = [1;1];
% eta0_AP_02 = [1;1];
% eta0_PL_02 = [1];
% eta0_DO_02 = [1;1];
% eta0_OT_02 = [1;1];
% eta0_CL_11 = [1;1];
% eta0_AP_11 = [1;1];
% eta0_PL_11 = [1];
% eta0_DO_11 = [1;1];
% eta0_OT_11 = [1;1];
% eta0_CL_20 = [1;1];
% eta0_AP_20 = [1;1];
% eta0_PL_20 = [1];
% eta0_DO_20 = [1;1];
% eta0_OT_20 = [1;1];
load('opt_initialize.mat')
for i =1:2
    i
    [etaCL,nllCL] = fmincon(@(eta)NLL_CL(eta,x),eta0_CL,[],[],[],[],[0;0],[pi;pi],[],options);
    [etaAP,nllAP] = fmincon(@(eta)NLL_AP(eta,x),eta0_AP,[],[],[],[],[0],[pi],[],options);
    [etaPL,nllPL] = fmincon(@(eta)NLL_DO(eta,x,0),eta0_PL,[],[],[],[],[1e-2],[pi],[],options);
    [etaDO,nllDO] = fmincon(@(eta)NLL_DO(eta,x,acos(0.5)),eta0_DO,[],[],[],[],[1e-2;1e-2],[pi;pi]/2,[],options);
    [etaOT,nllOT] = fmincon(@(eta)NLL_DO(eta,x,pi/2),eta0_OT,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    eta0_CL = etaCL
    eta0_AP = etaAP
    eta0_PL = etaPL
    eta0_DO = etaDO
    eta0_OT = etaOT
end
%%
for i =1:2
    i
    [etaCL,nllCL] = fmincon(@(eta)NLL_CL(eta,x_02),eta0_CL_02,[],[],[],[],[0;0],[pi;pi],[],options);
    [etaAP,nllAP] = fmincon(@(eta)NLL_AP(eta,x_02),eta0_AP_02,[],[],[],[],[0],[pi],[],options);
    [etaPL,nllPL] = fmincon(@(eta)NLL_DO(eta,x_02,0),eta0_PL_02,[],[],[],[],[1e-2],[pi],[],options);
    [etaDO,nllDO] = fmincon(@(eta)NLL_DO(eta,x_02,acos(0.5)),eta0_DO_02,[],[],[],[],[1e-2;1e-2],[pi;pi]/2,[],options);
    [etaOT,nllOT] = fmincon(@(eta)NLL_DO(eta,x_02,pi/2),eta0_OT_02,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    eta0_CL_02 = etaCL
    eta0_AP_02 = etaAP
    eta0_PL_02 = etaPL
    eta0_DO_02 = etaDO
    eta0_OT_02 = etaOT
end

for i =1:2
    i
    [etaCL,nllCL] = fmincon(@(eta)NLL_CL(eta,x_11),eta0_CL_11,[],[],[],[],[0;0],[pi;pi],[],options);
    [etaAP,nllAP] = fmincon(@(eta)NLL_AP(eta,x_11),eta0_AP_11,[],[],[],[],[0],[pi],[],options);
    [etaPL,nllPL] = fmincon(@(eta)NLL_DO(eta,x_11,0),eta0_PL_11,[],[],[],[],[1e-2],[pi],[],options);
    [etaDO,nllDO] = fmincon(@(eta)NLL_DO(eta,x_11,acos(0.5)),eta0_DO_11,[],[],[],[],[1e-2;1e-2],[pi;pi]/2,[],options);
    [etaOT,nllOT] = fmincon(@(eta)NLL_DO(eta,x_11,pi/2),eta0_OT_11,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    eta0_CL_11 = etaCL
    eta0_AP_11 = etaAP
    eta0_PL_11 = etaPL
    eta0_DO_11 = etaDO
    eta0_OT_11 = etaOT
end

for i =1:2
    i
    [etaCL,nllCL] = fmincon(@(eta)NLL_CL(eta,x_20),eta0_CL_20,[],[],[],[],[0;0],[pi;pi],[],options);
    [etaAP,nllAP] = fmincon(@(eta)NLL_AP(eta,x_20),eta0_AP_20,[],[],[],[],[0],[pi],[],options);
    [etaPL,nllPL] = fmincon(@(eta)NLL_DO(eta,x_20,0),eta0_PL_20,[],[],[],[],[1e-2],[pi],[],options);
    [etaDO,nllDO] = fmincon(@(eta)NLL_DO(eta,x_20,acos(0.5)),eta0_DO_20,[],[],[],[],[1e-2;1e-2],[pi;pi]/2,[],options);
    [etaOT,nllOT] = fmincon(@(eta)NLL_DO(eta,x_20,pi/2),eta0_OT_20,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    eta0_CL_20 = etaCL
    eta0_AP_20 = etaAP
    eta0_PL_20 = etaPL
    eta0_DO_20 = etaDO
    eta0_OT_20 = etaOT
end
%%
% seed_02 = zeros(1000, Ns_02);
% seed_11 = zeros(1000, Ns_11);
% seed_20 = zeros(1000, Ns_20);
% seed_combined = zeros(1000, Ns);
% for i = 1:1000
%     seed_02(i,:) = randi(Ns_02,[1,Ns_02]);
%     seed_11(i,:) = randi(Ns_11,[1,Ns_11]);
%     seed_20(i,:) = randi(Ns_20,[1,Ns_20]);
%     seed_combined(i,:) = randi(Ns,[1,Ns]);
% end
load('booststrap_seed.mat')
load('opt_initialize.mat')
%% Contralateral

BSNum = 800;
etaopt_CL_1  = zeros(2,BSNum);
etaopt_CL_2  = zeros(2,BSNum);
etaopt_CL_3  = zeros(2,BSNum);
etaopt_CL  = zeros(2,BSNum);
nllopt_CL = zeros(4,BSNum);
for i = 1:BSNum
    r_02 = seed_02(200+i,:);
    r_11 = seed_11(200+i,:);
    r_20 = seed_20(200+i,:);
    r = seed_combined(200+i,:);
    tic
    [etaopt_CL_1(:,i),nllopt_CL(1,i)] = fmincon(@(eta)NLL_CL(eta,x_02(:,r_02)),eta0_CL_02,[],[],[],[],[0;0],[pi;pi],[],options);

    [etaopt_CL_2(:,i),nllopt_CL(2,i)] = fmincon(@(eta)NLL_CL(eta,x_11(:,r_11)),eta0_CL_11,[],[],[],[],[0;0],[pi;pi],[],options);

    [etaopt_CL_3(:,i),nllopt_CL(3,i)] = fmincon(@(eta)NLL_CL(eta,x_20(:,r_20)),eta0_CL_20,[],[],[],[],[0;0],[pi;pi],[],options);

    [etaopt_CL(:,i),nllopt_CL(4,i)] = fmincon(@(eta)NLL_CL(eta,x(:,r)),eta0_CL,[],[],[],[],[0;0],[pi;pi],[],options);
    toc
    save("BS_CL_new2.mat","etaopt_CL","etaopt_CL_1","etaopt_CL_2","etaopt_CL_3","nllopt_CL");
end

%% Antipodal
BSNum = 800;
etaopt_AP_1  = zeros(2,BSNum);
etaopt_AP_2  = zeros(2,BSNum);
etaopt_AP_3  = zeros(2,BSNum);
etaopt_AP  = zeros(2,BSNum);
nllopt_AP = zeros(4,BSNum);
for i = 1:BSNum
    r_02 = seed_02(200+i,:);
    r_11 = seed_11(200+i,:);
    r_20 = seed_20(200+i,:);
    r = seed_combined(200+i,:);
    tic
    [etaopt_AP_1(:,i),nllopt_AP(1,i)] = fmincon(@(eta)NLL_AP(eta,x_02(:,r_02)),eta0_AP_02,[],[],[],[],[0],[pi],[],options);

    [etaopt_AP_2(:,i),nllopt_AP(2,i)] = fmincon(@(eta)NLL_AP(eta,x_11(:,r_11)),eta0_AP_11,[],[],[],[],[0],[pi],[],options);

    [etaopt_AP_3(:,i),nllopt_AP(3,i)] = fmincon(@(eta)NLL_AP(eta,x_20(:,r_20)),eta0_AP_20,[],[],[],[],[0],[pi],[],options);

    [etaopt_AP(:,i),nllopt_AP(4,i)] = fmincon(@(eta)NLL_AP(eta,x(:,r)),eta0_AP,[],[],[],[],[0],[pi],[],options);
    toc
    save("BS_AP_new2.mat","etaopt_AP","etaopt_AP_1","etaopt_AP_2","etaopt_AP_3","nllopt_AP");
end
%% Orthogonal
BSNum = 800;
% etaopt_OT_1  = zeros(2,BSNum);
% etaopt_OT_2  = zeros(2,BSNum);
% etaopt_OT_3  = zeros(2,BSNum);
% etaopt_OT  = zeros(2,BSNum);
% nllopt_OT = zeros(4,BSNum);
for i = 746:BSNum
    r_02 = seed_02(200+i,:);
    r_11 = seed_11(200+i,:);
    r_20 = seed_20(200+i,:);
    r = seed_combined(200+i,:);
    tic
    [etaopt_OT_1(:,i),nllopt_OT(1,i)] = fmincon(@(eta)NLL_DO(eta,x_02(:,r_02),pi/2),eta0_OT_02,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);

    [etaopt_OT_2(:,i),nllopt_OT(2,i)] = fmincon(@(eta)NLL_DO(eta,x_11(:,r_11),pi/2),eta0_OT_11,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);

    [etaopt_OT_3(:,i),nllopt_OT(3,i)] = fmincon(@(eta)NLL_DO(eta,x_20(:,r_20),pi/2),eta0_OT_20,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);

    [etaopt_OT(:,i),nllopt_OT(4,i)] = fmincon(@(eta)NLL_DO(eta,x(:,r),pi/2),eta0_OT,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    toc
    save("BS_OT_new2.mat","etaopt_OT","etaopt_OT_1","etaopt_OT_2","etaopt_OT_3","nllopt_OT");
end
% Parallel
BSNum = 800;
chi = 0;
etaopt_PL_1  = zeros(3,BSNum);
etaopt_PL_2  = zeros(3,BSNum);
etaopt_PL_3  = zeros(3,BSNum);
etaopt_PL  = zeros(3,BSNum);
nllopt_PL = zeros(4,BSNum);
for i = 1:BSNum
    r_02 = seed_02(200+i,:);
    r_11 = seed_11(200+i,:);
    r_20 = seed_20(200+i,:);
    r = seed_combined(200+i,:);
    tic
    [etaopt_PL_1(:,i),nllopt_PL(1,i)] = fmincon(@(eta)NLL_DO(eta,x_02(:,r_02),chi),eta0_PL_02,[],[],[],[],[1e-2],[pi],[],options);

    [etaopt_PL_2(:,i),nllopt_PL(2,i)] = fmincon(@(eta)NLL_DO(eta,x_11(:,r_11),chi),eta0_PL_11,[],[],[],[],[1e-2],[pi],[],options);

    [etaopt_PL_3(:,i),nllopt_PL(3,i)] = fmincon(@(eta)NLL_DO(eta,x_20(:,r_20),chi),eta0_PL_20,[],[],[],[],[1e-2],[pi],[],options);

    [etaopt_PL(:,i),nllopt_PL(4,i)] = fmincon(@(eta)NLL_DO(eta,x(:,r),chi),eta0_PL,[],[],[],[],[1e-2],[pi],[],options);
    toc
    save("BS_PL_new2.mat","etaopt_PL","etaopt_PL_1","etaopt_PL_2","etaopt_PL_3","nllopt_PL");
end

% distance-optimal with u/v=0.5
BSNum = 800;
chi = acos(0.5);
etaopt_DO_1  = zeros(2,BSNum);
etaopt_DO_2  = zeros(2,BSNum);
etaopt_DO_3  = zeros(2,BSNum);
etaopt_DO  = zeros(2,BSNum);
nllopt_DO = zeros(4,BSNum);
for i = 1:BSNum
    r_02 = seed_02(200+i,:);
    r_11 = seed_11(200+i,:);
    r_20 = seed_20(200+i,:);
    r = seed_combined(200+i,:);
    tic
    [etaopt_DO_1(:,i),nllopt_DO(1,i)] = fmincon(@(eta)NLL_DO(eta,x_02(:,r_02),chi),eta0_DO_02,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);

    [etaopt_DO_2(:,i),nllopt_DO(2,i)] = fmincon(@(eta)NLL_DO(eta,x_11(:,r_11),chi),eta0_DO_11,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);

    [etaopt_DO_3(:,i),nllopt_DO(3,i)] = fmincon(@(eta)NLL_DO(eta,x_20(:,r_20),chi),eta0_DO_20,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);

    [etaopt_DO(:,i),nllopt_DO(4,i)] = fmincon(@(eta)NLL_DO(eta,x(:,r),chi),eta0_DO,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    toc
    save("BS_DO_new2.mat","etaopt_DO","etaopt_DO_1","etaopt_DO_2","etaopt_DO_3","nllopt_DO");
end
%% varying chi
x_20_psi = x_20;
x_20_psi(2,:) = wrapToPi(x_20(2,:)+x_20(1,:)+pi);
load('opt_initialize.mat')
chi_list = linspace(pi/36,pi/2,18);
nll_storage = zeros(1,length(chi_list));
eta_storage = zeros(3,length(chi_list));
eta0_DO(1) = eta0_PL;
eta0_DO(2) = eta0_PL;
eta0_DO(3) = eta0_PL;
save('pre_results_20_psiphi.mat',"eta_storage","nll_storage")
for i =1:length(chi_list)
    i
    [etaDO,nllDO] = fmincon(@(eta)NLL_DO(eta,x_20_psi,chi_list(i)),eta0_DO,[],[],[],[],[2e-2;2e-2;2e-2],[pi;pi;pi],[],options);
    eta0_DO = etaDO
    nllDO
    nll_storage(i) = nllDO;
    eta_storage(:,i) = eta0_DO;
end
save('pre_results_20_psiphi.mat',"eta_storage","nll_storage")
%% for chi=0
load('opt_initialize.mat')

nll_bootstrap_all = zeros(200,1);
eta_bootstrap_all = zeros(200,1);
nll_bootstrap_02 = zeros(200,1);
eta_bootstrap_02 = zeros(200,1);
nll_bootstrap_11 = zeros(200,1);
eta_bootstrap_11 = zeros(200,1);
nll_bootstrap_20 = zeros(200,1);
eta_bootstrap_20 = zeros(200,1);
for k = 1:200
    r_all = randi(Ns,[1,Ns]);
    r_02 = randi(Ns_02,[1,Ns_02]);
    r_11 = randi(Ns_11,[1,Ns_11]);
    r_20 = randi(Ns_20,[1,Ns_20]);
    [etaPL_all,nllPL_all] = fmincon(@(eta)NLL_DO(eta,x(:,r_all),0),eta0_PL,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    nllPL
    nll_bootstrap_all(k) = nllPL_all;
    eta_bootstrap_all(k) = etaPL_all;
    [etaPL_02,nllPL_02] = fmincon(@(eta)NLL_DO(eta,x_02(:,r_all),0),eta0_PL,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    nllPL
    nll_bootstrap_02(k) = nllPL_02;
    eta_bootstrap_02(k) = etaPL_02;
    [etaPL_11,nllPL_11] = fmincon(@(eta)NLL_DO(eta,x_11(:,r_all),0),eta0_PL,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    nllPL
    nll_bootstrap_11(k) = nllPL_11;
    eta_bootstrap_11(k) = etaPL_11;
    [etaPL_20,nllPL_20] = fmincon(@(eta)NLL_DO(eta,x_20(:,r_all),0),eta0_PL,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
    nllPL
    nll_bootstrap_20(k) = nllPL_20;
    eta_bootstrap_20(k) = etaPL_20;

    save('bestchiresults_all_zero.mat', ...
        "eta_bootstrap_all","nll_bootstrap_all",...
        "eta_bootstrap_02","nll_bootstrap_02", ...
        "eta_bootstrap_11","nll_bootstrap_11", ...
        "eta_bootstrap_20","nll_bootstrap_20")
end
%% for chi>0
load('pre_results.mat')
chi_list = linspace(pi/36,pi/2,18);
nll_bootstrap = zeros(200,length(chi_list));
eta_bootstrap = zeros(200,2,length(chi_list));
for k = 1:200
    r_all = randi(Ns,[1,Ns]);
    for i =1:length(chi_list)
        eta0_DO = eta_storage(:,i);
        [etaDO,nllDO] = fmincon(@(eta)NLL_DO(eta,x(:,r_all),chi_list(i)),eta0_DO,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
        %     eta0_DO = etaDO
        nllDO
        nll_bootstrap(k,i) = nllDO;
        eta_bootstrap(k,:,i) = etaDO;
    end
    save('bestchiresults_all.mat',"eta_bootstrap","nll_bootstrap")
end
%%
load('pre_results_02.mat')
chi_list = linspace(pi/36,pi/2,18);
nll_bootstrap = zeros(200,length(chi_list));
eta_bootstrap = zeros(200,2,length(chi_list));
for k = 1:200
    r_02 = randi(Ns_02,[1,Ns_02]);
    for i =1:length(chi_list)
        eta0_DO = eta_storage(:,i);
        [etaDO,nllDO] = fmincon(@(eta)NLL_DO(eta,x_02(:,r_02),chi_list(i)),eta0_DO,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
        %     eta0_DO = etaDO
        nllDO
        nll_bootstrap(k,i) = nllDO;
        eta_bootstrap(k,:,i) = etaDO;
    end
    save('bestchiresults_02.mat',"eta_bootstrap","nll_bootstrap")
end
%%
load('pre_results_11.mat')
chi_list = linspace(pi/36,pi/2,18);
nll_bootstrap = zeros(200,length(chi_list));
eta_bootstrap = zeros(200,2,length(chi_list));
for k = 1:200
    r_11 = randi(Ns_11,[1,Ns_11]);
    for i =1:length(chi_list)
        eta0_DO = eta_storage(:,i);
        [etaDO,nllDO] = fmincon(@(eta)NLL_DO(eta,x_11(:,r_11),chi_list(i)),eta0_DO,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
        %     eta0_DO = etaDO
        nllDO
        nll_bootstrap(k,i) = nllDO;
        eta_bootstrap(k,:,i) = etaDO;
    end
    save('bestchiresults_11.mat',"eta_bootstrap","nll_bootstrap")
end
%%
load('pre_results_20.mat')
chi_list = linspace(pi/36,pi/2,18);
nll_bootstrap = zeros(200,length(chi_list));
eta_bootstrap = zeros(200,2,length(chi_list));
for k = 1:200
    r_20 = randi(Ns_20,[1,Ns_20]);
    for i =1:length(chi_list)
        eta0_DO = eta_storage(:,i);
        [etaDO,nllDO] = fmincon(@(eta)NLL_DO(eta,x_20(:,r_20),chi_list(i)),eta0_DO,[],[],[],[],[1e-2;1e-2],[pi;pi],[],options);
        %     eta0_DO = etaDO
        nllDO
        nll_bootstrap(k,i) = nllDO;
        eta_bootstrap(k,:,i) = etaDO;
    end
    save('bestchiresults_20.mat',"eta_bootstrap","nll_bootstrap")
end
%%
load('optNLLs.mat');
N = 120;

map_CL = zeros(N,N);
map_AP = zeros(N,N);


sig_theta = linspace(0.01, pi-0.01, N);
sig_phi = linspace(0.01, pi-0.01, N);

for i = 1:N
    tic
    sig_phi0 = sig_phi(i);
    for j = 1:N
        map_CL(i,j) = NLL_CL([sig_phi0,sig_theta(j)],x);
    end
    toc
end

for i = 1:N
    tic
    sig_phi0 = sig_phi(i);
    for j = 1:N
        map_AP(i,j) = NLL_AP([sig_phi0,sig_theta(j)],x);
    end
    toc
end

%%
load('optNLLs.mat');
N = 120;

map_PL = zeros(N,N);
map_DO = zeros(N,N);
map_OT = zeros(N,N);
map_PLsym = zeros(N,N);
map_DOsym = zeros(N,N);
map_OTsym = zeros(N,N);

sig_theta = linspace(0.01, pi-0.01, N);
sig_phi = linspace(0.01, pi-0.01, N);
sig_lam = linspace(0.01, pi-0.01, N);

sig_phi0 = etaPL(1);
for i = 1:N
    tic
    sig_lam0 = sig_lam(i);
    for j = 1:N
        map_PL(i,j) = NLL_DO([sig_phi0,sig_lam0,sig_theta(j)],x,0);
    end
    toc
end

sig_phi0 = etaDO(1);
for i = 1:N
    tic
    sig_lam0 = sig_lam(i);
    for j = 1:N
        map_DO(i,j) = NLL_DO([sig_phi0,sig_lam0,sig_theta(j)],x,acos(0.5));
    end
    toc
end

sig_phi0 = etaOT(1);
for i = 1:N
    tic
    sig_lam0 = sig_lam(i);
    for j = 1:N
        map_OT(i,j) = NLL_OT([sig_phi0,sig_lam0,sig_theta(j)],x);
    end
    toc
end

for i = 1:N
    tic
    sig_lam0 = sig_lam(i);
    for j = 1:N
        map_PLsym(i,j) = NLL_DO([sig_phi(j),sig_lam0,sig_theta(j)],x,0);
    end
    toc
end

for i = 1:N
    tic
    sig_lam0 = sig_lam(i);
    for j = 1:N
        map_DOsym(i,j) = NLL_DO([sig_phi(j),sig_lam0,sig_theta(j)],x,acos(0.5));
    end
    toc
end

for i = 1:N
    tic
    sig_lam0 = sig_lam(i);
    for j = 1:N
        map_OTsym(i,j) = NLL_OT([sig_phi(j),sig_lam0,sig_theta(j)],x);
    end
    toc
end


save mapDATA.mat map_PL map_DO map_OT map_PLsym map_DOsym map_OTsym
