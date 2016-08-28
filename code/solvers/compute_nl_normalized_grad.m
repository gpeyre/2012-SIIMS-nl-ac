function [Phi_fid Ener] =compute_nl_normalized_grad(HePhi, nx,ny,q,q1,w1,nl_dis,GauK)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Compute the shape gradient "-grad E" with E = E^{NL-N}
% 
% grad E = ||grad phi(x)|| * "- Phi_fid"
%         = ||grad phi(x)|| * [(Gx_in./Hx_in + Cx_in) - (Gx_out./Hx_out + Cx_out)]
% where
% Gx_in = int_{Om} K(x,y) dy  where K(x,y) = G_sigma(x,y)*d(px,py)
% Hx_in = int_{Om} G_sigma(x,y) dy
% Cx_in = int_{Om} [K(y,x)Hy_in - G_sigma(y,x)Gy_in]/(Hy_in)^2 dy
%
% Gx_out = int_{Om^c} K(x,y) dy  
% Hx_out = int_{Om^c} G_sigma(x,y) dy
% Cx_out = int_{Om^c} [K(y,x)Hy_out - G_sigma(y,x)Gy_out]/(Hy_out)^2 dy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('computing NL-Normalized model......');

G_in=zeros(nx,ny);G_out=zeros(nx,ny);
H_in=zeros(nx,ny);H_out=zeros(nx,ny);
C_in=zeros(nx,ny);C_out=zeros(nx,ny);

KD=zeros(q1,q1);  

%% compute Gx_in, Hx_in, Gx_out, Hx_out   
for x=1:nx           
     for y=1:ny
         selx = x-q:x+q;                
         sely = y-q:y+q;               
         selx(selx<1 | selx>nx) = [];
         sely(sely<1 | sely>ny) = [];
         
        %% K(x,y) = G_sigma(x,y)*d(px,py)
         KD(q1+1-length(selx):q1,q1+1-length(sely):q1)= nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1); 
        %% G_sigma(x,y)
         K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1,:);            

        %% Gx_in, Hx_in, Gx_out, Hx_out   
         G_in(x,y)= sum(sum((HePhi(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                
         H_in(x,y)= sum(sum((HePhi(selx,sely)).*K));                                            

         G_out(x,y)= sum(sum((1-HePhi(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); %% outside omega                
         H_out(x,y)= sum(sum((1-HePhi(selx,sely)).*K));                                       
     end           
end

%% compute Cx_in, Cx_out
for x=1:nx            
    for y=1:ny                
        selx = x-q:x+q;               
        sely = y-q:y+q;                
        selx(selx<1 | selx>nx) = [];               
        sely(sely<1 | sely>ny) = [];
        
       %% K(x,y) = G_sigma(x,y)*d(px,py)
        KD(q1+1-length(selx):q1,q1+1-length(sely):q1)= nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1); 
       %% G_sigma(x,y)
        K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1,:);     
        
       %% Cx_in, Cx_out  
        %%% Cx_in = int_{Om} [K(y,x)Hy_in - G_sigma(y,x)Gy_in]/(Hy_in)^2 dy) 
        A_in = (KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H_in(selx,sely)-K.*G_in(selx,sely))./(H_in(selx,sely).^2);              
        C_in(x,y)= sum(sum((HePhi(selx,sely)).*A_in)); 
        
        %%% Cx_out = int_{Om^c} [K(y,x)Hy_out - G_sigma(y,x)Gy_out]/(Hy_out)^2 dy        
        A_out = (KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H_out(selx,sely)-K.*G_out(selx,sely))./(H_out(selx,sely).^2);                
        C_out(x,y)= sum(sum((1-HePhi(selx,sely)).*A_out));              
    end    
end

%%% -grad E^{NL-U}
Phi_fid= - (G_in./H_in + C_in) + (G_out./H_out + C_out);
Ener = HePhi.*(G_in./H_in)+(1-HePhi).*(G_out./H_out);

Phi_fid=Phi_fid/(nx*ny*w1*w1);
Ener=Ener/(nx*ny*w1*w1);

