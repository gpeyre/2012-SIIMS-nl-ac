function [Phi_fid Ener]=compute_nl_unnormalized_grad(HePhi,nx,ny,q,q1,w1,nl_dis)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute the shape gradient of "- grad E" with E = E^{NL-U}
%  grad E = ||grad phi(x)|| * "Phi_fid"
%         = ||grad phi(x)|| * 
%            (int_{Om} K(x,y) + K(y,x) dy - int_{Om^c} K(x,y) + K(y,x) dy)
% where K(x,y) = G_sigma(x,y)*d(px,py).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('computing NL-Un-Normalized model......');

A1=zeros(nx,ny); 
A2=zeros(nx,ny);

KD=zeros(q1,q1);
for x=1:nx
        for y=1:ny
            selx = x-q:x+q;
            sely = y-q:y+q;
            
            selx(selx<1 | selx>nx) = [];
            sely(sely<1 | sely>ny) = [];

           %% K(x,y) = G_sigma(x,y)*d(px,py)
            KD(q1+1-length(selx):q1,q1+1-length(sely):q1)= nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1);
            
           %% int_{Om} K(x,y) dy
            A1(x,y)=sum(sum(HePhi(selx,sely).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
           %% int_{Om^c} K(x,y) dy
            A2(x,y)=sum(sum((1-HePhi(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                                 
        end      
end
Phi_fid = -2*A1+2*A2; %% due to the symmetry K(x,y) = K(y,x)
% Phi_fid = 2*A1-2*A2; %% due to the symmetry K(x,y) = K(y,x)
Ener=HePhi.*A1+(1-HePhi).*A2;

Phi_fid=Phi_fid/(nx*ny*q1*q1*w1*w1);
Ener=Ener/(nx*ny*q1*q1*w1*w1);
