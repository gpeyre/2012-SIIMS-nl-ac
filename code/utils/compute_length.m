function norm_gradHePhi=compute_length(phi, Heavieps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the length of curve C={x: phi(x) = 0}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gradHePhi = grad(Heaviside_eps(phi,Heavieps));
norm_gradHePhi= sqrt(sum(gradHePhi.^2,3));
norm_gradHePhi=norm_gradHePhi.*(norm_gradHePhi<1)+norm_gradHePhi.*(norm_gradHePhi==1)+1.*(norm_gradHePhi>1);