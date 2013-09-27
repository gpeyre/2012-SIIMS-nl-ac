function [phi Energy] = sol_gac(phi, M, Edge, options)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Geodesic Active Contours (GAC)
% 
% E^{GAC} = int_{C} g(x) dx where g(x) is an edge function
% Compute LSF evolution with balloon force term "eta*Edge*|nabla phi|"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('solving GAC model......');

[nx,ny]=size(M); n=nx;

Heavieps = options.Heavieps;
niter = options.niter;
redis_num = options.redis_num;
dt = options.dt;
eta = options.eta;
tol = options.tol;

options.order=2;

figure;
for iter=1:niter
    iter
     
     %% update phi
     gD = grad(phi,options);
     d1 = max(eps, sqrt(sum(gD.^2,3)) );
     d = sqrt(sum(gD.^2,3));
     g = gD ./ repmat( d1, [1 1 2] );
     
     G = d.*div( repmat(Edge,[1 1 2]).*g,options )/n+eta*Edge.*d;
     phi = phi + dt*G;
    
    %% re-distancing 
    if mod(iter,redis_num)==0
        phi = perform_redistancing(phi);
    end
    
    %% compute Energy functional  
    norm_gradHePhi=compute_length(phi, Heavieps); 
    Energy(iter) = sum(sum(Edge.*norm_gradHePhi))/n;  
    
    %% stopping
    if iter>1
        if sqrt((Energy(iter)-Energy(iter-1)).^2)/sqrt(Energy(iter).^2) <tol
            break;
        end
    end       

   %% plot 
    plot_levelset(phi,0,M);title(['Iteration # ',int2str(iter)]);
    mov(iter) = getframe;  
    
end
