function [phi Energy] = sol_cv_wasser_color(phi, M, MC, P, w1, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Color image => (1) Brightness + (2) Chromaticity
% 
% * Decompose a color image M to (1) B = |M| & (2) C = M/|M|
% * Apply 1d-Wasserstein metric to (1) Brightness component (B)
% * Apply Multichannel Chan-Vese model to (2) Chromaticity (C)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('solving CV-Wasser-color......');

[nx ny L]=size(M); n=nx;

Heavieps = options.Heavieps;
niter = options.niter;
redis_num = options.redis_num;
dt = options.dt;
alpha= options.alpha;
beta= options.beta;
gamma = options.gamma;
tol = options.tol;

options.order=2;
c1=zeros(1,L); c2=zeros(1,L);

figure; 
for iter=1:niter
    iter
     
    %% Brightness: Patch (P) + Wasserstein distance
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
       
    for x=1:nx
        for y=1:ny    
            Fid_tmp = sum((P(x,y,:)- P1).^2,3); 
            Fid_tmp2 = sum((P(x,y,:)- P2).^2,3);            
            Fid(x,y)= (Fid_tmp2 - Fid_tmp)/(nx*ny*w1*w1);
            
            tmp1 = sum((P(x,y,:)- P1).^2.*(HePhipat(x,y,:)),3); 
            tmp2= sum((P(x,y,:)- P2).^2.*(1-HePhipat(x,y,:)),3); 
            Ener(x,y) = (tmp1+tmp2)/(nx*ny*w1*w1); 
        end
    end

   %% Chromaticity component: multichannel-CV 
   ain=sum(sum(HePhi));
   aout=sum(sum(1-HePhi));
   for l=1:3
        c1(l)=sum(sum((HePhi).*MC(:,:,l)))/ain;
        c2(l)=sum(sum((1-HePhi).*MC(:,:,l)))/aout;   
    end
    Fid2=zeros(n,n);
    for l=1:3
         Fid2=Fid2+ (MC(:,:,l)-c1(l)).^2 -(MC(:,:,l)-c2(l)).^2;
    end
    Fid2 = Fid2/(nx*ny);
         
    %% Update "phi"                       
    options.order=2;
    gD = grad(phi,options);       
    d = max(eps, sqrt(sum(gD.^2,3)) );
    g = gD ./ repmat( d, [1 1 2] );
     
    G = d.*(gamma*div(g,options)/n-beta*Fid2+alpha*Fid);  
    phi= phi+ dt* G;  
                     
    %% re-distancing 
    if mod(iter,redis_num)==0  
       phi = perform_redistancing(phi);
    end
    
    %% compute Energy functional
    Ener1= sum(Ener(:));
    A1=zeros(nx,ny); A2=zeros(nx,ny);
    for l=1:L
        A1=A1+(MC(:,:,l)-c1(l)).^2;
        A2=A2+(MC(:,:,l)-c2(l)).^2;
    end   
    Ener2 = sum(sum(A1.*HePhi+A2.*(1-HePhi)))/(nx*ny);
    norm_gradHePhi=compute_length(phi, Heavieps);
    Ener3 = sum(sum(norm_gradHePhi))/n; 
    Energy(iter) = alpha*Ener1+beta*Ener2+gamma*Ener3;
    
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
