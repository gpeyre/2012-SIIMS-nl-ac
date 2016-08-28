function nl_dis = compute_motion_dist(M, M2, w1, w, q1, q,  GauK, metric, options)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Motion profiles Pr(x) & Motion distance d(Pr(x), Pr(y)) 
% metric = 'motion-1d' (gray-scaled), 'motion-3d' (color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx ny L] = size(M);

sigm = options.sigm;
ndx = options.ndx;
ndy = options.ndy;
w2= options.w2;

%% compute K(x,y) = G_sigma(x,y)*d(Px,Py)
switch metric  
    case 'motion-1d'        
        nl_dis=compute_motion_profile_distance_gray(M,M2, nx, ny, w,w1,w2, q, q1, sigm, ndx, ndy, GauK);
        
    case 'motion-3d'
        nl_dis=compute_motion_profile_distance_color(M,M2, nx, ny, w,w1, w2, q,q1, sigm, ndx, ndy, GauK);
end