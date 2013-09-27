function [phi phi2 Energy]=sol_multiphase_cv_MI(M, phi, phi2, options)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiphase MI-Chan-Vese energy (Chan-Vese Multiphase model)
%
% Chan-Vese Multiphase model with "m=2" level set function 
% (up to N=2^m regions can be segmented)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('solving CV-multiphase-MI......');

[nx ny]=size(M); n=nx;

Heavieps = options.Heavieps;
niter = options.niter;
redis_num = options.redis_num;
dt = options.dt;
lambda = options.lambda;
gamma = options.gamma;
tol = options.tol;

figure;
for itr=1:niter
    itr
    
    %% update means "c1, c2" 
    HePhi = Heaviside_eps(phi,Heavieps);
    HePhi2 = Heaviside_eps(phi2,Heavieps);
    
    a1=sum(sum(HePhi.*HePhi2));
    a2=sum(sum(HePhi.*(1-HePhi2)));    
    a3=sum(sum((1-HePhi).*HePhi2));
    a4=sum(sum((1-HePhi).*(1-HePhi2)));
    
    c1=sum(sum(HePhi.*HePhi2.*M))/a1;
    c2=sum(sum(HePhi.*(1-HePhi2).*M))/a2;   
    c3=sum(sum((1-HePhi).*HePhi2.*M))/a3;  
    c4=sum(sum((1-HePhi).*(1-HePhi2).*M))/a4;  
    
    %% grad E^{MI-CV} 
    Fid1= ((M-c1).^2-(M-c3).^2).*(HePhi2)+ ((M-c2).^2-(M-c4).^2).*(1-HePhi2);
    Fid2= ((M-c1).^2-(M-c2).^2).*(HePhi)+ ((M-c3).^2-(M-c4).^2).*(1-HePhi);
     
    %% - grad L(phi)
    options.order=2;
    gD = grad(phi,options);
    d1 = max(eps, sqrt(sum(gD.^2,3)) );
    d = sqrt(sum(gD.^2,3)) ;
    g = gD ./ repmat( d1, [1 1 2] );
    
    %% update phi1, phi2 : phi = phi + dt (-grad L(phi) - grad E^{MI-CV})
    G1 = d.*(gamma*div(g,options)/n-lambda*(Fid1)/(nx*ny)); 
    G2 = d.*(gamma*div(g,options)/n-lambda*(Fid2)/(nx*ny));
    
    phi = phi + dt*G1;
    phi2 = phi2 + dt*G2;
    
    %% Re-distancing 
    if mod(itr,redis_num)==0
        phi = perform_redistancing(phi);
        phi2 = perform_redistancing(phi2);
    end
    
    %% Energy functional    
    E1 = (M-c1).^2.*HePhi.*HePhi2+(M-c2).^2.*HePhi.*(1-HePhi2);
    E2 = (M-c3).^2.*(1-HePhi).*HePhi2+(M-c4).^2.*(1-HePhi).*(1-HePhi2);
    Ener1 = sum(E1(:)+E2(:))/(nx*ny);
    norm_gradHePhi=compute_length(phi, Heavieps);
    Ener2 = sum(sum(norm_gradHePhi))/n;   
    Energy(itr) = lambda*Ener1+gamma*Ener2;
    
    %% stopping
    if itr>1
        if sqrt((Energy(itr)-Energy(itr-1)).^2)/sqrt(Energy(itr).^2) <tol
            break;
        end
    end       
    
    %% Plot 
    imagesc(M);colormap(gray);axis image off;hold on;contour(phi,[0 0], 'r','linewidth',3); 
    hold on;contour(phi2,[0 0], 'y','linewidth',3); 
    mov(itr) = getframe;   

end

