function [phi Energy] = sol_other_methods(phi, M, M0, method, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GAC, CV, IAC (pixel based)
% Gray-valued, Color, Gabor transforms
% method = 'gac', 'cv' (1d, multichannel), 'iac' (1d, multichannel)
%
%% Extensions of CV, LBF model incorporated with Patches
% Gray-valued image/textures (L2-1d / Wasserstein-1d)
% method = 'cv-patch', 'lbf-patch'
%
%% 'cv-wasser-color':
% Color images/textures: wasserstein + multichannel chan-vese model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx ny L] = size(M);

%%
% color image M => Brightness (MB) + Chromaticity (MC) components 
switch method 
case 'cv-wasser-color'    
    MB=sqrt(M(:,:,1).^2+M(:,:,2).^2+M(:,:,3).^2);
    MC=zeros(nx,ny,3);
    for i=1:3
    MC(:,:,i)=M(:,:,i)./sqrt(M(:,:,1).^2+M(:,:,2).^2+M(:,:,3).^2);
    end
end

%% 
% extract patches for CV-patch, LBF-patch, cv-wasser-color
switch method 
case {'cv-patch','lbf-patch', 'cv-wasser-color'}
    
    w = options.w; w1 = 2*w+1;
    q = options.q; q1 = 2*q+1;
    sigma = options.sigma;
    metric = options.metric;
    
    switch method 
    case {'cv-patch','lbf-patch'}            
        P = extract_patch(M,nx,ny,w1,w);
    case 'cv-wasser-color'
        P = extract_patch(MB,nx,ny,w1,w);    
    end
    
    switch metric
        case 'wasserstein-1d'
        fprintf('extracting sorted patches......');
        for i=1:nx
            for j=1:ny
                  Pat(1:w1,1:w1)=P(i,j,1:w1,1:w1);
                  Pat_tmp=Pat;
                  Pat_sort(1:w1*w1,1) = sort(Pat_tmp(:));
                  Pnew(i,j,1:w1*w1) = Pat_sort(1:w1*w1,1);
            end
        end      
        P = Pnew;
    end
end

%% 
% windowing function G_sigma for LBF-patch
switch method 
    case 'lbf-patch'
    Gaussi=fspecial('gaussian',[q1 q1],sigma); 
    GauK=Gaussi/max(Gaussi(:));  
    GauKp=zeros(q1,q1,w1*w1);        
    for x=1:q1             
        for y=1:q1 
            GauKp(x,y,1:w1*w1)=GauK(x,y); 
        end    
    end
end

%% 
switch method
    case 'gac'
        Edge = options.Edge;
        [phi Energy] = sol_gac(phi, M, Edge, options);  
    case {'cv','iac'}
        [phi Energy] = sol_cv_iac(phi, M, M0, method, options);    
    case 'cv-patch'
        [phi Energy] = sol_cv_patch(phi,M,P,w1,options);
    case 'lbf-patch'
        [phi Energy]= sol_lbf_patch(phi,M,P,w1, q, q1, GauKp, GauK, options); 
    case 'cv-wasser-color'
       [phi Energy] = sol_cv_wasser_color(phi, M, MC, P, w1, options);
end

