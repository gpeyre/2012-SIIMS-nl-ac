function nl_dis = compute_l2_wasser_dist(M, w1, w, q1, q, GauK, metric, options)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L2 distance (weighted)
% metric = 'l2-1d' (gray-scaled), 'l2-3d' (color), 'l2-gabor' (Gabor transforms)
%
%% Wasserstein distance
% metric = 'wasserstein-1d' (gray-scaled), 'wasserstein-3d' (Sliced Wasserstein, color)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx ny L] = size(M);

%% compute K(x,y) = G_sigma(x,y)*d(Px,Py)
switch metric
    case 'l2-1d'
       gau_a = options.gau_a;
       nl_dis=compute_l2_dist_gray(M,nx,ny,w1,w,gau_a,q,q1,GauK);       
    case 'l2-3d'
       gau_a = options.gau_a;
       nl_dis=compute_l2_dist_color(M,nx,ny,w1,w,gau_a,q,q1,GauK);       
    case 'wasserstein-1d'    
       nl_dis=compute_wasser_dist_gray(M,nx,ny,w1,w,q,q1,GauK);      
    case 'wasserstein-3d'
       nl_dis=compute_wasser_dist_color(M,nx,ny,w1,w,q,q1,GauK);        
    case 'l2-gabor'
       gau_a = options.gau_a; 
       Mg=compute_gabor_transforms(M,nx,ny,options); %% gabor transformed features
       nl_dis=compute_l2_dist_gabor(Mg,nx,ny,w1,w,gau_a,q,q1,GauK); 
end
