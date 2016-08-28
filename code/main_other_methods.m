%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GAC, CV, IAC (pixel based)
% Gray-valued, Color, Gabor transforms
% @ method = 'gac', 'cv' (1d, multichannel), 'iac' (1d, multichannel)
% @ compute_edge = 'edge-1d', 'edge-3d', 'edge-gabor' (only for method ='gac','iac')
% @ parameters:: eta ('gac'), gamma ('cv'), mu ('iac')
% @ Gabor transforms (for example='l2-gabor'):: 
%   sigma_start, [freq_start, freq_end], ntheta, nsigma, nfeq 
%
%% Extensions of CV, LBF model incorporated with Patches
% Gray-valued image/textures (L2-1d / Wasserstein-1d)
% @ method = 'cv-patch', 'lbf-patch'
% @ parameters:: 
% * w (half size of patch "p", w1=2*w+1)
% * q, sigma (half size of G_sigma(x,y), q1=2*q+1)
% * gamma (regularization parameter) 
%
%% 'cv-wasser-color':
% Color images/textures: wasserstein + multichannel chan-vese model
% @ parameters:: alpha (wasserstein term), beta (Multi-CV term), gamma (regularization) 
%
%% Multiphase models: MR (Samson's et al), MI (Chan-Vese model)
% method = 'multiphase-1d-mr' (2 or 3 LSFs), 'multiphase-1d-mi' (2 LSFs)
% @ parameters:: gamma (regularization), mu (only for MR, repulsive term)
%
% Common Parameters (for the LSF evolution)::
% * niter (# of iterations), * redis_num (redistancing step for LSF)
% * dt (step size for the gradient descent method)
%
% ### October, 2011, Miyoun Jung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% example = 'multiphase-1d-mr'; % 2 or 3 LSFs
% example = 'multiphase-1d-mi'; % 2 LSFs

[M M2 phi phi2 phi3]= load_image_and_initial_lsf(example);
[nx ny L]=size(M); n=nx; M0 = M;

figure;imagesc(M);colormap(gray);axis image off;  hold on;contour(phi,[0 0], 'r','linewidth',4); 
% hold on;contour(phi2,[0 0], 'c','linewidth',4); 
% hold on;contour(phi3,[0 0], 'y','linewidth',4); 
figure;mesh(phi); 
% hold on; mesh(phi2); 
% hold on; mesh(phi3);

%% choose a method 
method = 'gac';
method = 'cv';
% method = 'iac';
% method = 'cv-patch';
% method = 'lbf-patch';
% method = 'cv-wasser-color';
% method = 'samson-multiphase-mr';
% method = 'cv-multiphase-mi';

%% choose parameters
options.Heavieps = 1.0000e-010; 
options.Heavieps = 1; 
options.lambda = 1; % fixed

options.niter = 500;
options.redis_num = 100;
options.tol = 1.0000e-010; 

switch method    
case {'cv-patch', 'lbf-patch', 'cv-wasser-color'}
    options.w = 5;                       % half size of patch
    options.q = 15;                      % half size of G_sigma
    options.sigma = 1000;                % G_sigma (= 1000/n)
    options.metric = 'l2-1d';            % choose metric
    options.metric = 'wasserstein-1d';
    w1 = 2*options.w+1;
    q1 = 2*options.q+1;
end

switch method
case {'cv','iac'}
    options.dt = 2000*n;       % CV, IAC
    options.gamma = 0.005/n;   % CV
    options.mu = 0.3/n;        % IAC
case {'gac'}
    options.dt =100;         
    options.eta = -0.1/n;      
case 'cv-patch'
    options.dt = 10*(n*w1*w1);
    options.gamma = 5/(n*w1*w1); 
case 'lbf-patch'
    options.dt = 0.03*(n*q1*q1*w1*w1);
    options.gamma = 2000/(n*q1*q1*w1*w1); 
case 'cv-wasser-color'
    options.dt = 2*(n*w1*w1);
    options.gamma = 0.2/(n*w1*w1); 
    options.alpha = 10;
    options.beta = 1/(w1*w1);
case {'samson-multiphase-mr', 'cv-multiphase-mi'}
    options.dt = 300*n;
    options.gamma = 0.1/n;
    options.mu = 0.1;  % only for MR method
end

%% Gabor transforms (Only for L2-Gabor)
switch example
case 'l2-gabor'
    options.sigma_start = 4;
    options.freq_start = 2; options.freq_end = 6;
    options.ntheta = 1; 
    options.nsigma = 2; 
    options.nfreq = 2;
    Mg = compute_gabor_transforms(M,nx,ny,options); 
    M=Mg; [nx ny L]=size(M); 
end

%% compute Edge function
switch method
case {'gac', 'iac'}    
    compute_edge='edge-1d'; 
    % compute_edge='edge-3d';
    % compute_edge='edge-gabor';
    
    options.alpha = 1;
    options.epsilon = 0.1;
    options.gau_sig = 2; % (= 0.5/n = 2/4*n), gau_sig/4 = usual sigma
    Edge = edge_detector(M, compute_edge, options);
    options.Edge = Edge;
end

%% solve other methods 
switch method
case {'gac','cv','iac','cv-patch','lbf-patch','cv-wasser-color'}
    [phi Energy] = sol_other_methods(phi, M, M0, method ,options);
case 'samson-multiphase-mr'
    choose_num ='2';
    [phi phi2 phi3 Energy]=sol_multiphase_samson_MR(M, phi, phi2, phi3, options, choose_num);
case 'cv-multiphase-mi'
    [phi phi2 Energy]=sol_multiphase_cv_MI(M, phi, phi2, options);
end

return;

figure;imagesc(M0);colormap(gray);axis image off;hold on;
contour(phi,[0 0],'r','linewidth',4);
