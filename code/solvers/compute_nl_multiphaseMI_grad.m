%% choose a function: depending on # LSF = 2 or 3

% function [Phi_fid Phi_fid2 Phi_fid3 Ener]=compute_nl_multiphaseMI_grad(nx,ny,q,q1,w1,nl_dis,GauK,HePhi, HePhi2,HePhi3, nl_model, choose_num)  
function [Phi_fid Phi_fid2 Ener]=compute_nl_multiphaseMI_grad(nx,ny,q,q1,w1,nl_dis,GauK,HePhi, HePhi2, nl_model,choose_num)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the shape gradient "-grad E" of NL-multiphase-MI-model
% 
% 1. choose_nl_model: "un-normalized" (NL-U) or "normalized" (NL-N)
% 2. choose number of level set functions: "2" or "3"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=@(u,v)u.*v+(1-u).*(1-v);

KD=zeros(q1,q1);
Phi_fid=zeros(nx,ny);
Phi_fid2=zeros(nx,ny);
Phi_fid3=zeros(nx,ny);

switch nl_model 
    
   %% un-normalized (NL-U)
    case 'unnormalized'        
              
          for x=1:nx
            for y=1:ny
                selx = x-q:x+q;
                sely = y-q:y+q;

                selx(selx<1 | selx>nx) = [];
                sely(sely<1 | sely>ny) = [];

                %%% K(x,y) = G_sigma(x,y) d(px,py)
                KD(q1+1-length(selx):q1,q1+1-length(sely):q1)= nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1);

                U = repmat(HePhi(x,y), [length(selx) length(sely)]); % H(phi)
                V = HePhi(selx,sely); % H(-phi)
                U2 = repmat(HePhi2(x,y), [length(selx) length(sely)]); % H(phi2)
                V2 = HePhi2(selx,sely); % H(-phi2)
                  
                switch choose_num   
                    case '2'  
                   %% phi1
                     A1=sum(sum(HePhi(selx,sely).*g(U2,V2).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                     A2=sum(sum((1-HePhi(selx,sely)).*g(U2,V2).*(KD(q1+1-length(selx):q1,q1+1-length(sely):q1)))); 
                   %% phi2
                     B1=sum(sum(HePhi2(selx,sely).*g(U,V).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                     B2=sum(sum((1-HePhi2(selx,sely)).*g(U,V).*(KD(q1+1-length(selx):q1,q1+1-length(sely):q1)))); 
                    
                     Phi_fid(x,y) = -2*A1+2*A2; 
                     Phi_fid2(x,y) = -2*B1+2*B2;
                     Ener(x,y) = sum(sum(g(U,V).*g(U2,V2).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                 
                     case '3' 
                     U3 = repmat(HePhi3(x,y), [length(selx) length(sely)]); % H(phi3)
                     V3 = HePhi3(selx,sely); % H(-phi3)
                     
                   %% phi1
                     A1=sum(sum(HePhi(selx,sely).*g(U2,V2).*g(U3,V3).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                     A2=sum(sum((1-HePhi(selx,sely)).*g(U2,V2).*g(U3,V3).*(KD(q1+1-length(selx):q1,q1+1-length(sely):q1)))); 
                   %% phi2
                     B1=sum(sum(HePhi2(selx,sely).*g(U,V).*g(U3,V3).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                     B2=sum(sum((1-HePhi2(selx,sely)).*g(U,V).*g(U3,V3).*(KD(q1+1-length(selx):q1,q1+1-length(sely):q1)))); 
                   %% phi3
                     C1=sum(sum(HePhi3(selx,sely).*g(U,V).*g(U2,V2).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); 
                     C2=sum(sum((1-HePhi3(selx,sely)).*g(U,V).*g(U2,V2).*(KD(q1+1-length(selx):q1,q1+1-length(sely):q1)))); 

                     Phi_fid(x,y) = -2*A1+2*A2; 
                     Phi_fid2(x,y) = -2*B1+2*B2;
                     Phi_fid3(x,y) = -2*C1+2*C2;
                     Ener(x,y) = sum(sum(g(U,V).*g(U2,V2).*g(U3,V3).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                    
                end
            end
          end          
          Phi_fid=Phi_fid/(nx*ny*q1*q1*w1*w1); 
          Phi_fid2=Phi_fid2/(nx*ny*q1*q1*w1*w1); 
          switch choose_num 
          case '3' 
            Phi_fid3 = Phi_fid3/(nx*ny*q1*q1*w1*w1); 
          end
          Ener=Ener/(nx*ny*q1*q1*w1*w1);
                             
    case 'normalized' 
        
    %% normalized (NL-N) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % CASE. Two level set functions phi1, phi2::
     % For phi1, 
     % (1) g(u,v) = u*v + (1-u)*(1-v) with u=H(phi2(x)), v=H(phi2(y))
     % (2) K^{1}(x,y)=K(x,y) * g(u,v), G^{1}_sigma(x,y)=G_sigma(x,y)*g(u,v)
     %     where K(x,y) = G_sigma(x,y)*d(px,py)
     %
     % Gx1_in = int_{Om1} K^{1}(x,y) dy  
     % Hx1_in = int_{Om1} G^{1}_sigma(x,y) * g(u,v) dy
     % Cx1_in = int_{Om1} [K^{1}(y,x)Hy_in - G_sigma(y,x)Gy_in]/(Hy_in)^2 dy
     %
     % Gx1_out = int_{Om1^c} K^{1}(x,y) dy  
     % Hx1_out = int_{Om1^c} K^{1}_sigma(x,y) dy
     % Cx1_out = int_{Om1^c} [K^{1}(y,x) Hy_out - K^{1}_sigma(y,x)Gy_out]/(Hy_out)^2 dy
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      %% compute Gx_in, Hx_in, Gx_out, Hx_out        
       for x=1:nx               
         for y=1:ny
             selx = x-q:x+q;                
             sely = y-q:y+q;               
             selx(selx<1 | selx>nx) = [];
             sely(sely<1 | sely>ny) = [];         

             %%% K(x,y) = G_sigma(x,y)d(px,py)
             KD(q1+1-length(selx):q1,q1+1-length(sely):q1)= nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1); 
             %%% G_sigma(x,y)
             K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1,:);          

             U = repmat(HePhi(x,y), [length(selx) length(sely)]);                
             V = HePhi(selx,sely);

             U2 = repmat(HePhi2(x,y), [length(selx) length(sely)]);               
             V2 = HePhi2(selx,sely);

             switch choose_num
             case '2'
               %% phi1   
                 G1_in(x,y)= sum(sum((HePhi(selx,sely)).*g(U2,V2).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                
                 H1_in(x,y)= sum(sum((HePhi(selx,sely)).*g(U2,V2).*K));                                            

                 G1_out(x,y)= sum(sum((1-HePhi(selx,sely)).*g(U2,V2).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); %% outside omega                
                 H1_out(x,y)= sum(sum((1-HePhi(selx,sely)).*g(U2,V2).*K));  
                  
               %% phi2
                 G2_in(x,y)= sum(sum((HePhi2(selx,sely)).*g(U,V).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                
                 H2_in(x,y)= sum(sum((HePhi2(selx,sely)).*g(U,V).*K));                                            

                 G2_out(x,y)= sum(sum((1-HePhi2(selx,sely)).*g(U,V).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); %% outside omega                
                 H2_out(x,y)= sum(sum((1-HePhi2(selx,sely)).*g(U,V).*K));          

             case'3'
                 U3 = repmat(HePhi3(x,y), [length(selx) length(sely)]);                   
                 V3 = HePhi3(selx,sely);

               %% phi1
                 G1_in(x,y)= sum(sum((HePhi(selx,sely)).*g(U2,V2).*g(U3,V3).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                
                 H1_in(x,y)= sum(sum((HePhi(selx,sely)).*g(U2,V2).*g(U3,V3).*K));                                            

                 G1_out(x,y)= sum(sum((1-HePhi(selx,sely)).*g(U2,V2).*g(U3,V3).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); %% outside omega                
                 H1_out(x,y)= sum(sum((1-HePhi(selx,sely)).*g(U2,V2).*g(U3,V3).*K));  

               %% phi2
                 G2_in(x,y)= sum(sum((HePhi2(selx,sely)).*g(U,V).*g(U3,V3).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                
                 H2_in(x,y)= sum(sum((HePhi2(selx,sely)).*g(U,V).*g(U3,V3).*K));                                            

                 G2_out(x,y)= sum(sum((1-HePhi2(selx,sely)).*g(U,V).*g(U3,V3).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); %% outside omega                
                 H2_out(x,y)= sum(sum((1-HePhi2(selx,sely)).*g(U,V).*g(U3,V3).*K));  

               %% phi3
                 G3_in(x,y)= sum(sum((HePhi3(selx,sely)).*g(U,V).*g(U2,V2).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1)));                
                 H3_in(x,y)= sum(sum((HePhi3(selx,sely)).*g(U,V).*g(U2,V2).*K));                                            

                 G3_out(x,y)= sum(sum((1-HePhi3(selx,sely)).*g(U,V).*g(U2,V2).*KD(q1+1-length(selx):q1,q1+1-length(sely):q1))); %% outside omega                
                 H3_out(x,y)= sum(sum((1-HePhi3(selx,sely)).*g(U,V).*g(U2,V2).*K));  
             end
         end         
       end
        
     %% compute Cx_in, Cx_out 
      for x=1:nx                         
        for y=1:ny                
            selx = x-q:x+q;               
            sely = y-q:y+q;                
            selx(selx<1 | selx>nx) = [];               
            sely(sely<1 | sely>ny) = [];

            KD(q1+1-length(selx):q1,q1+1-length(sely):q1)= nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1);                
            K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1,:);  
            
             U = repmat(HePhi(x,y), [length(selx) length(sely)]);                
             V = HePhi(selx,sely);

             U2 = repmat(HePhi2(x,y), [length(selx) length(sely)]);               
             V2 = HePhi2(selx,sely);

            switch choose_num
            case '2'
              %% phi1
                A1_in = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H1_in(selx,sely)-K.*G1_in(selx,sely)).*g(U2,V2))./(H1_in(selx,sely).^2);              
                C1_in(x,y)= sum(sum((HePhi(selx,sely)).*A1_in)); 

                A1_out = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H1_out(selx,sely)-K.*G1_out(selx,sely)).*g(U2,V2))./(H1_out(selx,sely).^2);                
                C1_out(x,y)= sum(sum((1-HePhi(selx,sely)).*A1_out));     

              %% phi2  
                A2_in = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H2_in(selx,sely)-K.*G2_in(selx,sely)).*g(U,V))./(H2_in(selx,sely).^2);              
                C2_in(x,y)= sum(sum((HePhi2(selx,sely)).*A2_in)); 

                A2_out = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H2_out(selx,sely)-K.*G2_out(selx,sely)).*g(U,V))./(H2_out(selx,sely).^2);                
                C2_out(x,y)= sum(sum((1-HePhi2(selx,sely)).*A2_out));     

            
            case '3'            
                U3 = repmat(HePhi3(x,y), [length(selx) length(sely)]);                   
                V3 = HePhi3(selx,sely);

              %% phi1   
                A1_in = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H1_in(selx,sely)-K.*G1_in(selx,sely)).*g(U2,V2).*g(U3,V3))./(H1_in(selx,sely).^2);              
                C1_in(x,y)= sum(sum((HePhi(selx,sely)).*A1_in)); 

                A1_out = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H1_out(selx,sely)-K.*G1_out(selx,sely)).*g(U2,V2).*g(U3,V3))./(H1_out(selx,sely).^2);                
                C1_out(x,y)= sum(sum((1-HePhi(selx,sely)).*A1_out));     

              %% phi2
                A2_in = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H2_in(selx,sely)-K.*G2_in(selx,sely)).*g(U,V).*g(U3,V3))./(H2_in(selx,sely).^2);              
                C2_in(x,y)= sum(sum((HePhi2(selx,sely)).*A2_in)); 

                A2_out = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H2_out(selx,sely)-K.*G2_out(selx,sely)).*g(U,V).*g(U3,V3))./(H2_out(selx,sely).^2);                
                C2_out(x,y)= sum(sum((1-HePhi2(selx,sely)).*A2_out));   

              %% phi3  
                A3_in = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H3_in(selx,sely)-K.*G3_in(selx,sely)).*g(U,V).*g(U2,V2))./(H3_in(selx,sely).^2);              
                C3_in(x,y)= sum(sum((HePhi3(selx,sely)).*A3_in)); 

                A3_out = ((KD(q1+1-length(selx):q1,q1+1-length(sely):q1).*H3_out(selx,sely)-K.*G3_out(selx,sely)).*g(U,V).*g(U2,V2))./(H3_out(selx,sely).^2);                
                C3_out(x,y)= sum(sum((1-HePhi3(selx,sely)).*A3_out));    
            end
        end  % end of for x
      end % end of for y
              
      Phi_fid= - (G1_in./H1_in + C1_in) + (G1_out./H1_out + C1_out);        
      Phi_fid2= - (G2_in./H2_in + C2_in) + (G2_out./H2_out + C2_out);
      Ener = HePhi.*(G1_in./H1_in)+(1-HePhi).*(G1_out./H1_out)+HePhi2.*(G2_in./H2_in)+(1-HePhi2).*(G2_out./H2_out);
      switch choose_num 
          case '3' 
            Phi_fid3 = - (G3_in./H3_in + C3_in) + (G3_out./H3_out + C3_out);
            Ener = HePhi.*(G1_in./H1_in)+(1-HePhi).*(G1_out./H1_out)+HePhi2.*(G2_in./H2_in)+(1-HePhi2).*(G2_out./H2_out)...
                +HePhi3.*(G3_in./H3_in)+(1-HePhi3).*(G3_out./H3_out);
            Phi_fid3 = Phi_fid3/(nx*ny*w1*w1);
      end      
      Phi_fid=Phi_fid/(nx*ny*w1*w1);
      Phi_fid2=Phi_fid2/(nx*ny*w1*w1);      
      Ener=Ener/(nx*ny*w1*w1);

end