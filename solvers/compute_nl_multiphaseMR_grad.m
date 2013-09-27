%% choose function: depending on # LSF = 2 or 3

% function [Phi_fid Phi_fid2 Phi_fid3 Ener]=compute_nl_multiphaseMR_grad(nx,ny,q,q1,w1,nl_dis,GauK,HePhi, HePhi2,HePhi3, nl_model, choose_num) 
function [Phi_fid Phi_fid2 Ener]=compute_nl_multiphaseMR_grad(nx,ny,q,q1,w1,nl_dis,GauK,HePhi, HePhi2,nl_model, choose_num)   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the shape gradient "-grad E" of NL-multiphase-MR-model
% 
% 1. choose_nl_model: "un-normalized" (NL-U) or "normalized" (NL-N)
% 2. choose number of level set functions: "2" or "3"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KD=zeros(q1,q1); 
Phi_fid=zeros(nx,ny);
Phi_fid2=zeros(nx,ny);
Phi_fid3=zeros(nx,ny);

switch nl_model
       
    %% unnormalized (NL-U)
    % (int_{Om} K(x,y) + K(y,x) dy - int_{Om^c} K(x,y) + K(y,x) dy)
    case 'unnormalized'        
        
        for x=1:nx
            for y=1:ny
                selx = x-q:x+q;
                sely = y-q:y+q;

                selx(selx<1 | selx>nx) = [];
                sely(sely<1 | sely>ny) = [];

                %%% K(x,y)=G_sigma(x,y)d(px,py) 
                KD(q1+1-length(selx):q1,q1+1-length(sely):q1)= nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1);
        
              %% phi1
                A1=sum(sum(HePhi(selx,sely).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                A2=sum(sum((1-HePhi(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
            
              %% phi2
                B1=sum(sum(HePhi2(selx,sely).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                B2=sum(sum((1-HePhi2(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                                                
                Phi_fid(x,y) = -2*A1+2*A2;
                Phi_fid2(x,y) = -2*B1+2*B2;
                
                Ener1(x,y)=A1; Ener2(x,y)=A2;
                Ener3(x,y)=B1; Ener4(x,y)=B2; 
                
                switch choose_num
                    case '3'
                  %% phi3   
                    C1=sum(sum(HePhi3(selx,sely).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                    C2=sum(sum((1-HePhi3(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                    Phi_fid3(x,y) = -2*C1+2*C2;
                    Ener5(x,y)=C1; Ener6(x,y)=C2;  
                end
            end 
        end
        Phi_fid=Phi_fid/(nx*ny*q1*q1*w1*w1); 
        Phi_fid2=Phi_fid2/(nx*ny*q1*q1*w1*w1); 
        Ener=(HePhi.*Ener1+(1-HePhi).*Ener2)+(HePhi2.*Ener3+(1-HePhi2).*Ener4);
        
        switch choose_num
            case '3'
            Phi_fid3=Phi_fid3/(nx*ny*q1*q1*w1*w1);
            Ener=(HePhi.*Ener1+(1-HePhi).*Ener2)+(HePhi2.*Ener3+(1-HePhi2).*Ener4)+(HePhi3.*Ener5+(1-HePhi3).*Ener6);
        end       
        Ener=Ener/(nx*ny*q1*q1*w1*w1);
        
      
    %% normalized (NL-N)
     % Remind:::
     % Gx_in = int_{Om} K(x,y) dy  where K(x,y) = G_sigma(x,y)*d(px,py)
     % Hx_in = int_{Om} G_sigma(x,y) dy
     % Cx_in = int_{Om} [K(y,x)Hy_in - G_sigma(y,x)Gy_in]/(Hy_in)^2 dy
     %
     % Gx_out = int_{Om^c} K(x,y) dy  
     % Hx_out = int_{Om^c} G_sigma(x,y) dy
     % Cx_out = int_{Om^c} [K(y,x)Hy_out - G_sigma(y,x)Gy_out]/(Hy_out)^2 dy
    case 'normalized'
           
            G_in=zeros(nx,ny);G_out=zeros(nx,ny);H_in=zeros(nx,ny);H_out=zeros(nx,ny);
            G_in2=zeros(nx,ny);G_out2=zeros(nx,ny);H_in2=zeros(nx,ny);H_out2=zeros(nx,ny);
            G_in3=zeros(nx,ny);G_out3=zeros(nx,ny);H_in3=zeros(nx,ny);H_out3=zeros(nx,ny);
            C_in=zeros(nx,ny);C_out=zeros(nx,ny);C_in2=zeros(nx,ny);C_out2=zeros(nx,ny);
            C_in3=zeros(nx,ny);C_out3=zeros(nx,ny);
            
            KD=zeros(q1,q1);   
           
          %% compute Gx_in, Hx_in, Gx_out, Hx_out for phi1, phi2, (or phi3)   
            for x=1:nx           
                 for y=1:ny
                     selx = x-q:x+q;                
                     sely = y-q:y+q;               
                     selx(selx<1 | selx>nx) = [];
                     sely(sely<1 | sely>ny) = [];

                     %%% K(x,y)=G_sigma(x,y)d(px,py) 
                     KD(q1+1-length(selx):q1,q1+1-length(sely):q1)= nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1); 
                     %%% G_sigma(x,y)
                     K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1,:);           

                  %% phi1   
                     G_in(x,y)= sum(sum((HePhi(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                
                     H_in(x,y)= sum(sum((HePhi(selx,sely)).*K));                                            

                     G_out(x,y)= sum(sum((1-HePhi(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); %% outside omega                
                     H_out(x,y)= sum(sum((1-HePhi(selx,sely)).*K));       
                     
                  %% phi2    
                     G_in2(x,y)= sum(sum((HePhi2(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                
                     H_in2(x,y)= sum(sum((HePhi2(selx,sely)).*K));                                            

                     G_out2(x,y)= sum(sum((1-HePhi2(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); %% outside omega                
                     H_out2(x,y)= sum(sum((1-HePhi2(selx,sely)).*K));  
                                            
                     switch choose_num
                     case '3'
                  %% phi3   
                     G_in3(x,y)= sum(sum((HePhi3(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                
                     H_in3(x,y)= sum(sum((HePhi3(selx,sely)).*K));                                            

                     G_out3(x,y)= sum(sum((1-HePhi3(selx,sely)).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); %% outside omega                
                     H_out3(x,y)= sum(sum((1-HePhi3(selx,sely)).*K));  
                     end
                 end           
            end

           %% compute Cx_in, Cx_out for phi1, phi2, (or phi3)   
            for x=1:nx            
                for y=1:ny                
                    selx = x-q:x+q;               
                    sely = y-q:y+q;                
                    selx(selx<1 | selx>nx) = [];               
                    sely(sely<1 | sely>ny) = [];

                    KD(q1+1-length(selx):q1,q1+1-length(sely):q1)= nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1);                
                    K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1,:);     

                  %% phi1
                    A_in = (KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H_in(selx,sely)-K.*G_in(selx,sely))./(H_in(selx,sely).^2);              
                    C_in(x,y)= sum(sum((HePhi(selx,sely)).*A_in)); 

                    A_out = (KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H_out(selx,sely)-K.*G_out(selx,sely))./(H_out(selx,sely).^2);                
                    C_out(x,y)= sum(sum((1-HePhi(selx,sely)).*A_out));  
                    
                  %% phi2
                    A_in2 = (KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H_in2(selx,sely)-K.*G_in2(selx,sely))./(H_in2(selx,sely).^2);              
                    C_in2(x,y)= sum(sum((HePhi2(selx,sely)).*A_in2)); 

                    A_out2 = (KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H_out2(selx,sely)-K.*G_out2(selx,sely))./(H_out2(selx,sely).^2);                
                    C_out2(x,y)= sum(sum((1-HePhi2(selx,sely)).*A_out2));  
                    
                    switch choose_num
                    case '3'
                  %% phi3 
                    A_in3 = (KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H_in3(selx,sely)-K.*G_in3(selx,sely))./(H_in3(selx,sely).^2);              
                    C_in3(x,y)= sum(sum((HePhi3(selx,sely)).*A_in3)); 

                    A_out3 = (KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H_out3(selx,sely)-K.*G_out3(selx,sely))./(H_out3(selx,sely).^2);                
                    C_out3(x,y)= sum(sum((1-HePhi3(selx,sely)).*A_out3));  
                    end
                end    
            end
            Phi_fid= - (G_in./H_in + C_in) + (G_out./H_out + C_out); 
            Phi_fid2= - (G_in2./H_in2 + C_in2) + (G_out2./H_out2 + C_out2);           
            Ener1 = HePhi.*(G_in./H_in)+(1-HePhi).*(G_out./H_out);
            Ener2 = HePhi2.*(G_in2./H_in2)+(1-HePhi2).*(G_out2./H_out2);
            
            Phi_fid=Phi_fid/(nx*ny*w1*w1);
            Phi_fid2=Phi_fid2/(nx*ny*w1*w1);
            Ener=(Ener1+Ener2)/(nx*ny*w1*w1);
                    
            switch choose_num                   
            case '3'
            Phi_fid3= - (G_in3./H_in3 + C_in3) + (G_out3./H_out3 + C_out3);
            Ener3 = HePhi3.*(G_in3./H_in3)+(1-HePhi3).*(G_out3./H_out3);            
            Phi_fid3=Phi_fid3/(nx*ny*w1*w1);
            Ener=(Ener1+Ener2+Ener3)/(nx*ny*w1*w1);
            end            
end
                                

                                