function [phi Energy]= sol_lbf_patch(phi,M,P,w1, q, q1, GauK, Gaussi, options) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extension of LBF model using Patches (LBF-patch)
% instead of pixel values as in the original paper 
% "Implicit Active Contours Driven by Local Binary Fitting Energy",
%  C.Li et al, 2007
%
% E^{lbf-patch} = int int G_sigma(x,y) d(Py,P{1,x}) H(phi(y)) dy dx 
%                  +  int int G_sigma(x,y) d(Py,P{2,x}) (1-H(phi(y))) dy dx   
% where P{1,x} and P{2,x} are the averaged patches inside and outside the
% curve, which are updated at every iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('solving LBF-patch......');

[nx ny]=size(M); n=nx;

Heavieps = options.Heavieps;
niter = options.niter;
redis_num = options.redis_num;
dt = options.dt;
lambda = options.lambda;
gamma = options.gamma;
tol = options.tol;

P1=zeros(nx,ny,w1*w1);
P2=zeros(nx,ny,w1*w1);

figure;
for iter=1:niter
    iter 
    
    %% update patches p1, p2
    fprintf('updating p1, p2......');

    HePhi = Heaviside_eps(phi,Heavieps);     
    HePhipat = repmat(HePhi,[1 1 w1 w1]); 
 
    for x=1:nx
        for y=1:ny
            selx = x-q:x+q;
            sely = y-q:y+q;
            
            selx(selx<1 | selx>nx) = [];
            sely(sely<1 | sely>ny) = [];

            K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1,:);                
                                            
            Ptmp=sum(sum(HePhipat(selx,sely,:).*K.*P(selx,sely,:),1),2)./(sum(sum(K.*HePhipat(selx,sely,:),1),2)+eps);
            P1(x,y,1:w1*w1)=Ptmp(1,1,1:w1*w1);
            
            Ptmp2=sum(sum((1-HePhipat(selx,sely,:)).*K.*P(selx,sely,:),1),2)./(sum(sum(K.*(1-HePhipat(selx,sely,:)),1),2)+eps);
            P2(x,y,1:w1*w1)=Ptmp2(1,1,1:w1*w1);
        end       
    end
        
    %% - grad E^{LBF-patch}      
    fprintf('computing -grad E^{LBF-patch}......');
    for x=1:nx
        for y=1:ny
            selx = x-q:x+q;
            sely = y-q:y+q;
            
            selx(selx<1 | selx>nx) = [];
            sely(sely<1 | sely>ny) = [];
            
            K2 = Gaussi(q1+1-length(selx):q1,q1+1-length(sely):q1);  
  
            Dis1 = sum((repmat(P(x,y,:),[length(selx) length(sely) 1])- P1(selx,sely,:)).^2, 3); 
            Dis2 = sum((repmat(P(x,y,:),[length(selx) length(sely) 1])- P2(selx,sely,:)).^2, 3); 
            
            Fid1 = sum(sum(K2.*Dis1)); 
            Fid2 = sum(sum(K2.*Dis2)); 
            Phi_fid(x,y)= Fid2-Fid1;
            
            Dis1 = sum((repmat(P1(x,y,:),[length(selx) length(sely) 1])- P(selx,sely,:)).^2, 3); 
            Dis2 = sum((repmat(P2(x,y,:),[length(selx) length(sely) 1])- P(selx,sely,:)).^2, 3); 
            Ener(x,y) = sum(sum((1-HePhi(selx,sely)).*K2.*Dis2))+sum(sum((HePhi(selx,sely)).*K2.*Dis1)); 
        end      
    end
         
    %% update "phi"                       
    options.order=2;
    gD = grad(phi,options);       
    d = max(eps, sqrt(sum(gD.^2,3)) );
    g = gD ./ repmat( d, [1 1 2] );
     
    G = d.*(gamma*div(g,options)/n+lambda*Phi_fid/(nx*ny*q1*q1*w1*w1));  
    phi= phi+ dt* G;  
                 
    %% re-Distancing 
    if mod(iter,redis_num)==0  
       phi = perform_redistancing(phi);
    end
    
    %% compute Energy functional
    Ener1= sum(Ener(:))/(nx*ny*q1*q1*w1*w1);
    norm_gradHePhi=compute_length(phi, Heavieps);    
    Ener2= sum(sum(norm_gradHePhi))/n;
    Energy(iter) = lambda*Ener1+gamma*Ener2;
    
    %% stopping
    if iter>1
        if sqrt((Energy(iter)-Energy(iter-1)).^2)/sqrt(Energy(iter).^2) <tol
            break;
        end
    end       
    
    %% plot
    plot_levelset(phi,0,M);title(['Iteration #',int2str(iter)]);
    mov(iter) = getframe;
           
end
