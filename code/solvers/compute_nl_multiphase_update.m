%% choose a function: depending on # LSF = 2 or 3

% function [phi phi2 phi3]=compute_nl_multiphase_update(phi,phi2,phi3,Phi_fid,Phi_fid2,Phi_fid3, nx,ny, multiphase_model, choose_num, dt, lambda, gamma, mu)
function [phi phi2]=compute_nl_multiphase_update(phi,phi2,Phi_fid,Phi_fid2, nx,ny, multiphase_model, choose_num, dt, lambda, gamma, mu)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Update LSFs of NL-multiphase-models (MR or MI):
% 
% Update phi, phi2, phi3 with given - grad E^{MR-NL,MI-NL}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.order=2; n=nx;
Heavieps =1.0000e-010;

HePhi = Heaviside_eps(phi,Heavieps);
HePhi2 = Heaviside_eps(phi2,Heavieps);
switch choose_num
case '3'
HePhi3 = Heaviside_eps(phi3,Heavieps);
end

switch multiphase_model
    
   %% NL-Multiphase-Repulsive model
    case 'MR'
        
    gD = grad(phi,options);       
    d = max(eps, sqrt(sum(gD.^2,3)) );
    g = gD ./ repmat( d, [1 1 2] );
    
    gD2 = grad(phi2,options);       
    d2 = max(eps, sqrt(sum(gD2.^2,3)) );
    g2 = gD2 ./ repmat( d2, [1 1 2] );
        
    switch choose_num
        case '2'
        F1=(1-(2*(HePhi2)))/(nx*ny);   
        F2=(1-(2*(HePhi)))/(nx*ny); 

        case'3'
        F1=(1-(2*(HePhi2+HePhi3)))/(nx*ny);
        F2=(1-(2*(HePhi+HePhi3)))/(nx*ny);
        F3=(1-(2*(HePhi+HePhi2)))/(nx*ny);
    end
     
    G = d.*(gamma*div(g,options)/n+lambda* Phi_fid+mu*F1);  
    phi= phi+ dt* G;  
    
    G2 = d2.*(gamma*div(g2,options)/n+lambda* Phi_fid2+mu*F2);  
    phi2= phi2+ dt* G2;
    
    switch choose_num
        case '3'
        gD3 = grad(phi3,options);       
        d3 = max(eps, sqrt(sum(gD3.^2,3)) );
        g3 = gD3 ./ repmat( d3, [1 1 2] );

        G3 = d3.*(gamma*div(g3,options)/n+lambda* Phi_fid3+mu*F3);  
        phi3= phi3+ dt* G3; 
    end
        
    %% NL-Multiphase-Intersection-model
    case 'MI'
    
    gD = grad(phi,options);       
    d = max(eps, sqrt(sum(gD.^2,3)) );
    g = gD ./ repmat( d, [1 1 2] );

    gD2 = grad(phi2,options);       
    d2 = max(eps, sqrt(sum(gD2.^2,3)) );
    g2 = gD2 ./ repmat( d2, [1 1 2] );
    
    G = d.*(gamma*div(g,options)/n+lambda* Phi_fid);  
    phi= phi+ dt* G;  
    
    G2 = d2.*(gamma*div(g2,options)/n+lambda* Phi_fid2);  
    phi2= phi2+ dt* G2; 
    
    switch choose_num
        case '3'
        gD3 = grad(phi3,options);       
        d3 = max(eps, sqrt(sum(gD3.^2,3)) );
        g3 = gD3 ./ repmat( d3, [1 1 2] );
        G3 = d3.*(gamma*div(g3,options)/n+lambda* Phi_fid3);  
        phi3= phi3+ dt* G3;
    end
        
end