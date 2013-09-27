%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non-local Active Contours
% ***Two NL models: "un-normalized" (NL-U) / " normalized" (NL-N)
%
% (1) L2-distance: pixel values (intensity, color, Gabor coefficient)
% (2) Wasserstein distance: texture 
% (3) Motion distance using motion profiles
%
% *** Multiphase method: MR (Repulsive) / MI (Intersection) method
% 
%% choose method:
% @ choose_LSF_num = '1', '2', '3'
% @ metric = 'l2-1d','l2-3d','l2-gabor','wasserstein-1d',
% 'wasserstein-3d','motion-1d', 'motion-3d'
% @ multiphase_model = 'MR', 'MI' (when Choose_LSF_num=2 or 3)
% @ nl_model = 'unnormalized', 'normalized'
%
%% choose (common) parameters 
% @ w: half size of patch "p" (or "pi" for motion distance)
% @ (q. sigma): half size of G_sigma(x,y)
% @ niter = total iteration # for LSF evolution
% @ redis_num = redistancing step for LSF 
% @ dt = step size for LSF evolution
% @ gamma = regularization parameter
%
%% Additional parameters
% @ L2-distance: gau_a = weights in the distance
% @ L2-Gabor: [ntheta, nsigma, nfreq] for Gabor transforms
%             sigma_start, [freq_start freq_end]
% @ Wasserstein-3d: choose the directions for projection 
%                   (here, fixed as [1,0,0],[0,1,0],[0,0,1])
% @ Motion-distance: ndx = half size of motion profile (=ndy)
%                    sigm = tuning parameter 
% @ Multiphase-MR method: mu = parameter for Repulsive term
%
% ### October,2011, Miyoun Jung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('toolbox/');
addpath('utils/');
addpath('solvers/');

close all
clear all

%% load an image and initial level set function
example ='l2-1d';
% example ='l2-gabor';
% example ='l2-3d';
% example ='wasserstein-1d';
% example ='wasserstein-3d';
% example = 'motion-1d';
% example = 'multiphase-1d-mr'; % 2 or 3 LSFs
% example = 'multiphase-1d-mi'; % 2 LSFs

[M M2 phi phi2 phi3]= load_image_and_initLSF(example);
[nx ny L]=size(M); n=nx;

figure;imagesc(M);colormap(gray);axis image off;  hold on;contour(phi,[0 0], 'r','linewidth',4); 
% hold on;contour(phi2,[0 0], 'c','linewidth',4); 
% hold on;contour(phi3,[0 0], 'y','linewidth',4); 
figure;mesh(phi); 
% hold on; mesh(phi2); 
% hold on; mesh(phi3); 

%% choose a metric & method 
choose_LSF_num ='1';
% choose_LSF_num ='2';
% choose_LSF_num ='3';

metric = 'l2-1d';
% metric = 'l2-3d';
% metric = 'l2-gabor';
% metric = 'wasserstein-1d';
% metric = 'wasserstein-3d';
% metric = 'motion-1d';
% metric = 'motion-3d';

multiphase_model = 'MR';
multiphase_model = 'MI';

nl_model = 'unnormalized';
% nl_model = 'normalized';

%% choose parameters 
w = 1;                         % half size of patch 
q = 15;                        % half size of G_sigma
sigma = 1000;                  % G_sigma (= 1000/n)

w1 = 2*w+1;                     % patch size w1*w1 
q1 = 2*q+1;                     % weight function support size q1*q1 
options.Heavieps = 1.0000e-010; % fixed
options.lambda = 1;             % fixed

options.nter = 100;               % number of iterations
options.redis_num = 10;           % re-distancing step
options.tol = 1.0000e-06;         % tolerance for stopping

switch nl_model  
case 'unnormalized'
    options.dt = 1*(n*q1*q1*w1*w1);
    options.gamma = 20/(n*q1*q1*w1*w1); 
    options.mu = 200/(q1*q1*w1*w1); % only for multiphase MR 
case 'normalized'
    options.dt = 100*(n*w1*w1);
    options.gamma = 0.05/(n*w1*w1); 
    options.mu = 200/(w1*w1); % only for multiphase MR 
end

switch metric
case {'l2-1d','l2-3d'}       
    options.gau_a = 0.5; %(= 0.5/n)
case 'l2-gabor'
    options.gau_a = 0.5; %(= 0.5/n)       
    options.sigma_start = 4;   
    options.freq_start = 2; options.freq_end = 6;        
    options.ntheta = 1;                  
    options.nsigma = 2;        
    options.nfreq = 2;
case {'motion-1d', 'motion-3d'}        
    options.sigm = 0.25; %(= 0.25/w1)
    options.ndx = 5; options.ndy = 5;       
    options.w2 = 0; %% fixed as 0 (local nbhd of pi(x))
    
    options.dt = 1000*n; %% normalized + motion
    options.gamma = 0.005/n; 
end

%% computing weights K(x,y) = G_sigma(x,y)d(px,py)
Gaussi=fspecial('gaussian',[q1 q1],sigma); 
GauK=Gaussi/max(Gaussi(:));  % windowing function G_sigma
switch metric
    case {'l2-1d','l2-3d','l2-gabor','wasserstein-1d','wasserstein-3d'}
        nl_dis = compute_l2_wasser_dist(M, w1, w, q1, q, GauK, metric, options);  
    case {'motion-1d', 'motion-3d'}
        nl_dis = compute_motion_dist(M, M2, w1, w, q1, q, GauK, metric, options); 
end

%% LSF evolution
switch choose_LSF_num 
    case '1'
        [phi Energy] = sol_nonlocal_seg(M, phi, q,q1,w1,nl_dis, GauK, nl_model, options);
    case '2'
        [phi phi2 Energy] = sol_nonlocal_seg_multiphase_lsf2(M, phi, phi2, q,q1,w1,nl_dis,GauK, multiphase_model, nl_model, options);
    case '3' 
        [phi phi2 phi3 Energy] = sol_nonlocal_seg_multiphase_lsf3(M,phi, phi2, phi3, q,q1,w1,nl_dis,GauK, multiphase_model, nl_model, options);  
end

% return;

figure;imagesc(M);colormap(gray);axis image off;  
hold on;contour(phi,[0 0], 'r','linewidth',4); 

figure;imagesc(M);colormap(gray);axis image off;  
hold on;contour(D0,[0 0], 'r','linewidth',4); 
