function [phi Energy] = sol_nonlocal_seg(M,phi,q,q1,w1,nl_dis, GauK, nl_model, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Non-local Active contours
% 1. nl_model = "un-normalized energy (NL-U)" or 
% 2. nl_model = "normalized energy (NL-N)" 
%
% Energy: E_total(phi) = E^{NL-U or NL-N} + gamma * Length(phi=0) 
% Evolution equation: phi = phi + dt * (- grad E(phi))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx ny L] = size(M); n=nx;

Heavieps = options.Heavieps;
lambda = options.lambda;
gamma = options.gamma;
niter = options.nter;
dt = options.dt;
redis_num = options.redis_num;
tol = options.tol;

figure; 
for iter=1:niter
    iter 

   %% compute "-grad E^(NL-U or NL-N)"    
    HePhi = Heaviside_eps(phi,Heavieps);   
    switch nl_model        
        case 'unnormalized'           
            [Phi_fid Ener]=compute_nl_unnormalized_grad(HePhi, nx,ny,q,q1,w1,nl_dis);        
        case 'normalized'           
            [Phi_fid Ener] =compute_nl_normalized_grad(HePhi, nx,ny,q,q1,w1,nl_dis,GauK);    
    end

   %% update "phi"    
    options.order=2;   
    gD = grad(phi,options);          
    d = max(eps, sqrt(sum(gD.^2,3)) );  
    g = gD ./ repmat( d, [1 1 2] );   
    G = d.*(gamma*div(g,options)/nx+lambda* Phi_fid); 

    phi= phi+ dt* G;  

    %% re-Distancing     
    if mod(iter,redis_num)==0  
           phi = perform_redistancing(phi);    
    end

   %% compute Energy functional
    Ener1(iter) = sum(Ener(:));
    norm_gradHePhi=compute_length(phi, Heavieps);
    Ener2(iter) = sum(sum(norm_gradHePhi))/n;    
    Energy(iter) =  lambda*Ener1(iter)+gamma*Ener2(iter);
    
%    %% exit
%     if iter>1
%         if abs(Energy(iter)-Energy(iter-1))/abs(Energy(iter)) < tol
%             break;
%         end
%     end
    
   %% plot
    plot_levelset(phi,0,M);title(['Iteration # ',int2str(iter)]);
    mov(iter) = getframe;
        
end
% movie(mov)
% fps=15;
% movie2avi(mov, 'ex.avi','fps',fps, 'compression','Cinepak');

