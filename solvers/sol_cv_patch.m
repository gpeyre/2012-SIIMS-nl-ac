function [phi Energy] = sol_cv_patch(phi,M,P,w1,options) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extension of Chan-Vese model using Patches 
%  instead of pixel values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('solving CV-patch......');

[nx ny]=size(M); n=nx;

Heavieps = options.Heavieps;
niter = options.niter;
redis_num = options.redis_num;
dt = options.dt;
lambda = options.lambda;
gamma = options.gamma;
tol = options.tol;

figure;
for iter=1:niter
    iter
     
    %% update patches p1, p2
     HePhi = Heaviside_eps(phi,Heavieps);      
     HePhipat = repmat(HePhi,[1 1 w1 w1]); 
    
     Ptmp=zeros(1,1,w1*w1);Ptmp2=zeros(1,1,w1*w1);
     Ptmp3=zeros(1,1,w1*w1);Ptmp4=zeros(1,1,w1*w1);
     for x=1:nx
        for y=1:ny                                               
            Ptmp=Ptmp + sum(sum(HePhipat(x,y,:).*P(x,y,:),1),2);
            Ptmp2=Ptmp2 +(sum(sum(HePhipat(x,y,:),1),2));
            
            Ptmp3=Ptmp3 +sum(sum((1-HePhipat(x,y,:)).*P(x,y,:),1),2);
            Ptmp4=Ptmp4 +(sum(sum((1-HePhipat(x,y,:)),1),2));
        end
     end
     P1=Ptmp./Ptmp2;
     P2=Ptmp3./Ptmp4;
       
     %% - grad E^{CV,patch}
     for x=1:nx
        for y=1:ny
            Fid1 = sum((P(x,y,:)- P1).^2,3); 
            Fid2 = sum((P(x,y,:)- P2).^2,3);            
            Phi_fid(x,y) = (Fid2-Fid1);
            
            tmp1 = sum((P(x,y,:)- P1).^2.*(HePhipat(x,y,:)),3); 
            tmp2= sum((P(x,y,:)- P2).^2.*(1-HePhipat(x,y,:)),3); 
            Ener(x,y) = (tmp1+tmp2); 
        end
     end
         
    %% update "phi"                       
    options.order=2;
    gD = grad(phi,options);       
    d = max(eps, sqrt(sum(gD.^2,3)) );
    g = gD ./ repmat( d, [1 1 2] );
     
    G = d.*(gamma*div(g,options)/n+lambda* Phi_fid/(nx*ny*w1*w1));  
    phi= phi+ dt* G;  
                    
    %% re-Distancing 
    if mod(iter,redis_num)==0  
       phi = perform_redistancing(phi);
    end
    
    %% compute Energy functional
    Ener1 = sum(Ener(:))/(nx*ny*w1*w1);
    norm_gradHePhi=compute_length(phi, Heavieps);    
    Ener2= sum(sum(norm_gradHePhi))/n;
    Energy(iter) = lambda*Ener1+gamma*Ener2;
    
   %% exit
    if iter>1
        if abs(Energy(iter)-Energy(iter-1))/abs(Energy(iter)) < tol
            break;
        end
    end
    
    %% plot
    plot_levelset(phi,0,M);title(['Iteration # ',int2str(iter)]);
    mov(iter) = getframe;
        
end
