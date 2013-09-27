function [phi phi2 Energy] = sol_nonlocal_seg_multiphase_lsf2(M, phi, phi2, q,q1,w1,nl_dis,GauK, multiphase_model, nl_model, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiphase MR or MI method with m=2 level set functions
% 
% nl-model: NL-U or NL-N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('solving multiphase model with 2 LSFs......');

[nx ny L] = size(M); n = nx;
choose_num = '2';

Heavieps = options.Heavieps;
lambda = options.lambda;
gamma = options.gamma;
mu = options.mu;
niter = options.nter;
dt = options.dt;
redis_num = options.redis_num;
tol = options.tol;

figure; 
for iter=1:niter
    iter 
    
    HePhi = Heaviside_eps(phi,Heavieps);
    HePhi2 = Heaviside_eps(phi2,Heavieps);
    
    %% - grad E^{MR, MI / NL-U, NL-N}
    switch multiphase_model 
        case 'MR'
           [Phi_fid Phi_fid2 Ener]=compute_nl_multiphaseMR_grad(nx,ny,q,q1,w1,nl_dis,GauK,HePhi, HePhi2,nl_model, choose_num);   
        case 'MI'
           [Phi_fid Phi_fid2 Ener]=compute_nl_multiphaseMI_grad(nx,ny,q,q1,w1,nl_dis,GauK,HePhi, HePhi2, nl_model,choose_num);    
    end
    
    %% Update "phi" 
    [phi phi2]=compute_nl_multiphase_update(phi,phi2,Phi_fid,Phi_fid2, nx,ny,multiphase_model, choose_num, dt, lambda, gamma, mu); 
        
    %% Re-Distancing 
    if mod(iter,redis_num)==0  
       phi = perform_redistancing(phi);
       phi2 = perform_redistancing(phi2);
    end
    
    %% compute Energy functional
    Ener1(iter) =  sum(Ener(:));
    switch multiphase_model
        case 'MR'
            F = (1-(HePhi+HePhi2)).^2;
            Ener2(iter) = sum(F(:))/(nx*ny);
    end   
    norm_gradHePhi=compute_length(phi, Heavieps);
    norm_gradHePhi2=compute_length(phi2, Heavieps);
    Ener3(iter) = sum(sum(norm_gradHePhi))/n+sum(sum(norm_gradHePhi2))/n;
    switch multiphase_model
        case 'MR'
        Energy(iter) = lambda*Ener1(iter)+mu*Ener2(iter)+gamma*Ener3(iter);
        case 'MI'
        Energy(iter) = lambda*Ener1(iter)+gamma*Ener3(iter);
    end
    
    %% stopping
    if iter>1
        if sqrt((Energy(iter)-Energy(iter-1)).^2)/sqrt(Energy(iter).^2) <tol
            break;
        end
    end       
    
    %% plot
    imagesc(M);colormap(gray);axis image off;hold on;contour(phi,[0 0], 'r','linewidth',3); 
    hold on;contour(phi2,[0 0], 'c','linewidth',3);
    mov(iter) = getframe;
        
end
