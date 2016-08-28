%% choose a function depending on # of LSFs
function [new_phi new_phi2 new_phi3 Energy]=sol_multiphase_samson_MR(M, phi, phi2, phi3, options, choose_num)
% function [phi phi2 Energy]=sol_multiphase_samson_MR(M, phi, phi2, options, choose_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiphase MR-Chan-Vese energy (Samson's Multiphase model)
%
% Samson's Multiphasemodel with "m=2" or "m=3" level set function 
% (N=m regions can be segmented)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fprintf('solving Samson-multiphase-MR......');

[nx ny]=size(M); n=nx;

Heavieps = options.Heavieps;
niter = options.niter;
redis_num = options.redis_num;
dt = options.dt;
lambda = options.lambda;
gamma = options.gamma;
mu = options.mu;
tol = options.tol;

figure;
for itr=1:niter
    itr
    
    %% update means c1, c2, c3
    HePhi = Heaviside_eps(phi,Heavieps);
    HePhi2 = Heaviside_eps(phi2,Heavieps);
    
    a1=sum(sum(HePhi));
    a2=sum(sum(HePhi2));    
    
    c1=sum(sum(HePhi.*M))/a1;
    c2=sum(sum(HePhi2.*M))/a2;   
    
    switch choose_num
    case '3'
    HePhi3 = Heaviside_eps(phi3,Heavieps);
    a3=sum(sum(HePhi3));
    c3=sum(sum(HePhi3.*M))/a3;  
    end
        
    %% - grad L(phi)
    options.order=2;
    gD = grad(phi,options);
    d1 = max(eps, sqrt(sum(gD.^2,3)) );
    d = sqrt(sum(gD.^2,3));
    g = gD ./ repmat( d1, [1 1 2] );
    
    gD2 = grad(phi2,options);       
    d2 = max(eps, sqrt(sum(gD2.^2,3)) );
    g2 = gD2 ./ repmat( d2, [1 1 2] );
    
    switch choose_num
    case '3'
    gD3 = grad(phi3,options);       
    d3 = max(eps, sqrt(sum(gD3.^2,3)));
    g3 = gD3 ./ repmat( d3, [1 1 2]);
    end
    
    %% Repulsive term: -grad (1- (Hphi+Hphi2+Hphi3))^2
    switch choose_num
    case '2'
        F = 2*(1-((HePhi+HePhi2)));
    case '3'
        F = 2*(1-((HePhi+HePhi2+HePhi3)));
    end

    %% update phi1, phi2, phi3
    G1 = d.*(gamma*div(g,options)/n-lambda*((M-c1).^2)/(nx*ny)+mu*F/(nx*ny)); 
    G2 = d2.*(gamma*div(g2,options)/n-lambda*((M-c2).^2)/(nx*ny)+mu*F/(nx*ny));
    phi = phi + dt*G1;
    phi2 = phi2 + dt*G2;
    
    switch choose_num
    case '2'
    phi3 = zeros(nx,ny);
    case '3'
    G3 = d3.*(gamma*div(g3,options)/n-lambda*((M-c3).^2)/(nx*ny)+2*mu*F/(nx*ny));
    phi3 = phi3 + dt*G3;
    end
           
    %% re-distancing 
    if mod(itr,redis_num)==0
        phi = perform_redistancing(phi);
        phi2 = perform_redistancing(phi2);       
        switch choose_num
        case '3'
        phi3 = perform_redistancing(phi3);
        end
    end
    
    %% compute Energy functional    
    switch choose_num
    case '2'
    Ener1 = sum(sum((M-c1).^2.*HePhi+(M-c2).^2.*HePhi2))/(nx*ny);
    norm_gradHePhi=compute_length(phi, Heavieps);
    norm_gradHePhi2=compute_length(phi2, Heavieps);
    Ener2 = sum(sum(norm_gradHePhi+norm_gradHePhi2))/n; 
    Ener3 = sum(sum((1- (HePhi+HePhi2)).^2))/(nx*ny);
    case '3'
    Ener1 = sum(sum((M-c1).^2.*HePhi+(M-c2).^2.*HePhi2+(M-c3).^2.*HePhi3))/(nx*ny);
    norm_gradHePhi=compute_length(phi, Heavieps);
    norm_gradHePhi2=compute_length(phi2, Heavieps);
    norm_gradHePhi3=compute_length(phi3, Heavieps);
    Ener2 = sum(sum(norm_gradHePhi+norm_gradHePhi2+norm_gradHePhi3))/n; 
    Ener3 = sum(sum((1- (HePhi+HePhi2+HePhi3)).^2))/(nx*ny);   
    end
    Energy(itr) =lambda*Ener1+gamma*Ener2+mu*Ener3;
    
    %% stopping
    if itr>1
        if sqrt((Energy(itr)-Energy(itr-1)).^2)/sqrt(Energy(itr).^2) <tol
            break;
        end
    end       
    
    %% plot 
    imagesc(M);colormap(gray);axis image off;hold on;contour(phi,[0 0], 'r','linewidth',3); 
    hold on;contour(phi2,[0 0], 'y','linewidth',3); 
    switch choose_num
    case '3'
    hold on;contour(phi3,[0 0], 'c','linewidth',3); 
    end
    mov(itr) = getframe;
    
end 
new_phi = phi; new_phi2 = phi2; new_phi3 = phi3;
