

Codes for the paper:

Nonlocal Active Contours (M. Jung, G. Peyré, L. D. Cohen), SIAM Journal on Imaging Sciences, vol. 5(3), pp. 1022-1054, 2012.
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% main_nonlocal_seg.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Non-local Active Contours
% ***Two Non-local models: "un-normalized" (NL-U) / " normalized" (NL-N)
%
% (1) L2-distance: pixel values (intensity, color, Gabor coefficient)
% (2) Wasserstein distance: texture 
% (3) Motion distance using motion profiles
%
% *** Multiphase methods: MR (Repulsive) / MI (Intersection) method
% 
%%% choose method:
% @ choose_LSF_num = '1', '2', '3'
% @ metric = 'l2-1d','l2-3d','l2-gabor','wasserstein-1d',
% 'wasserstein-3d','motion-1d', 'motion-3d'
% @ multiphase_model = 'MR', 'MI' (when Choose_LSF_num=2 or 3)
% @ nl_model = 'unnormalized', 'normalized'
%
%%% choose (common) parameters 
% @ w: half size of patch "p" (or "pi" for motion distance)
% @ (q. sigma): half size of G_sigma(x,y)
% @ niter = total iteration # for LSF evolution
% @ redis_num = redistancing step for LSF 
% @ dt = step size for LSF evolution
% @ gamma = regularization parameter
%
%%% Additional parameters
% @ L2-distance: gau_a = weights in the distance
% @ L2-Gabor: [ntheta, nsigma, nfreq] for Gabor transforms
%             sigma_start, [freq_start freq_end]
% @ Wasserstein-3d: choose the directions for projection 
%                   (here, fixed as [1,0,0],[0,1,0],[0,0,1])
% @ Motion-distance: ndx = half size of motion profile (=ndy)
%                    sigm = tuning parameter 
% @ Multiphase-MR method: mu = parameter for Repulsive term
%
%
% ### October,2011, Miyoun Jung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  main_other_methods.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GAC, CV, IAC (pixel based)
% Gray-valued, Color, Gabor transforms
% @ method = 'gac', 'cv' (1d, multichannel), 'iac' (1d, multichannel)
% @ compute_edge = 'edge-1d', 'edge-3d', 'edge-gabor' (only for method ='gac','iac')
% @ parameters:: eta ('gac'), gamma ('cv'), mu ('iac')
% @ Gabor transforms (for example='l2-gabor'):: 
%   sigma_start, [freq_start, freq_end], ntheta, nsigma, nfeq 
%
%%% Extensions of CV, LBF model incorporated with Patches
% Gray-valued image/textures (L2-1d / Wasserstein-1d)
% @ method = 'cv-patch', 'lbf-patch'
% @ parameters:: 
% * w (half size of patch "p", w1=2*w+1)
% * q, sigma (half size of G_sigma(x,y), q1=2*q+1)
% * gamma (regularization parameter) 
%
%%% 'cv-wasser-color':
% Color images/textures: wasserstein + multichannel chan-vese model
% @ parameters:: alpha (wasserstein term), beta (Multi-CV term), gamma (regularization) 
%
%%% Multiphase models: MR (Samson's et al), MI (Chan-Vese model)
% method = 'multiphase-1d-mr' (2 or 3 LSFs), 'multiphase-1d-mi' (2 LSFs)
% @ parameters:: gamma (regularization), mu (only for MR, repulsive term)
%
% Common Parameters (for the LSF evolution)::
% * niter (# of iterations), 
% * redis_num (redistancing step for LSF)
% * dt (step size for the gradient descent method)
%
%
% ### October, 2011, Miyoun Jung
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
