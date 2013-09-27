function [phi Energy] = sol_cv_iac(phi, M, M0, method, options)  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% "Chan-Vese" or "IAC" models
%
% method = 'cv', 'iac'
% IAC = Chan-Vese fidelity term + Geodesic length term (with edge detector)
% (gray-scaled/ color/ Gabor transforms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('solving CV or IAC models......');

[nx,ny,L]=size(M); n=nx;

Heavieps = options.Heavieps;
niter = options.niter;
redis_num = options.redis_num;
dt = options.dt;
lambda = options.lambda;
gamma = options.gamma;
mu= options.mu;
tol = options.tol;

options.order=2;
c1=zeros(1,L); c2=zeros(1,L);

figure;
for itr=1:niter
    itr
    
    %% update c1, c2
    HePhi = Heaviside_eps(phi,Heavieps);
    ain=sum(sum(HePhi));
    aout=sum(sum(1-HePhi));
    for l=1:L
        c1(l)=sum(sum((HePhi).*M(:,:,l)))/ain;
        c2(l)=sum(sum((1-HePhi).*M(:,:,l)))/aout;   
    end
      
    %% update phi 
    gD = grad(phi,options);
    d1 = max(eps, sqrt(sum(gD.^2,3)) );
    d = sqrt(sum(gD.^2,3)) ;
    g = gD ./ repmat( d1, [1 1 2] );    
    
    fid=zeros(nx,ny);
    for l=1:L
         fid=fid + lambda*(M(:,:,l)-c1(l)).^2 -lambda*(M(:,:,l)-c2(l)).^2;
    end
    
    switch method
        case 'cv'  %% Chan-Vese model
           G = d.*(gamma*div( g,options )/n-fid/(nx*ny));
   
        case 'iac' %% IAC model
           Edge = options.Edge;
           G = mu*d.*div(repmat(Edge,[1 1 2]).*g,options )/n-d.*(fid/(nx*ny));
    end
    phi = phi + dt*G;
    
    %% re-distancing 
    if mod(itr,redis_num)==0
        phi = perform_redistancing(phi);
    end
    
    %% energy functional    
    A1=zeros(nx,ny); A2=zeros(nx,ny);
    for l=1:L
        A1=A1+lambda*(M(:,:,l)-c1(l)).^2;
        A2=A2+lambda*(M(:,:,l)-c2(l)).^2;
    end   
    Ener1 = sum(sum(A1.*HePhi+A2.*(1-HePhi)))/(nx*ny);
    norm_gradHePhi=compute_length(phi, Heavieps);    
    Ener2= sum(sum(norm_gradHePhi))/n;
    Energy(itr) = lambda*Ener1+gamma*Ener2;
          
    %% exit
%     if itr>1
%         if abs(Energy(itr)-Energy(itr-1))/abs(Energy(itr)) < tol
%             break;
%         end
%     end

    %% plot 
    plot_levelset(phi,0,M0);title(['Iteration # ',int2str(itr)]);
    mov(itr) = getframe;    
end

