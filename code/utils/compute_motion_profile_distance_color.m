function nl_dis=compute_motion_profile_distance_color(M,M2, nx, ny, w,w1, w2, q,q1, sigm, ndx, ndy, GauK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Motion profiles Pr(x) and Motion distance d(Pr(x),Pr(y))
%% for color images
% 
% * Pr(x) is a motion profile at pixel x in the frame f^0: by comparing
%    - pi^{0}_x (a patch around x in first frame f^0=M) and
%    - pi^{1}_{x+de_i} (a patch around x+de_i in second frame f^1=M2) 
% 
% * d(Pr(x),Pr(y)) = 1 - sum_{i=1:tau*tau}( sqrt(Pr(x_i)) sqrt(Pr(y_i)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% extract patches from M (1st frame) and M2 (2nd frame)
P1r = extract_patch(M(:,:,1),nx,ny,w1,w);
P1g = extract_patch(M(:,:,2),nx,ny,w1,w);
P1b = extract_patch(M(:,:,3),nx,ny,w1,w);
P2r = extract_patch(M2(:,:,1),nx,ny,w1,w);
P2g = extract_patch(M2(:,:,2),nx,ny,w1,w);
P2b = extract_patch(M2(:,:,3),nx,ny,w1,w);

%% compute the motion profile
fprintf('computing Motion profile Pr......');

Pr=1+zeros(nx,ny,2*ndx+1,2*ndy+1);   
for x=1:nx
     for y=1:ny
         i=0; 
       %% range of "delta_i" in frame M2 (=[-ndx:ndx, -ndy:ndy])
        for dx=-ndx:ndx   
             i=i+1;j=0;            
             for dy=-ndy:ndy
                 j=j+1;                 
                selx = x-w2:x+w2; % fixed as w2=0
                sely = y-w2:y+w2;            
                selx(selx<1 | selx>nx | selx+dx<1 | selx+dx>nx) = [];
                sely(sely<1 | sely>ny | sely+dy<1 | sely+dy>ny) = [];
                
              %% patch distance |pi^0_x - pi^1_x+delta|^2
                D1 = (P1r(selx,sely,:) - P2r(selx+dx,sely+dy,:)).^2 ;                
                D2 = (P1g(selx,sely,:) - P2g(selx+dx,sely+dy,:)).^2;                 
                D3 = (P1b(selx,sely,:) - P2b(selx+dx,sely+dy,:)).^2 ;
                D = sum( (D1+D2+D3), 3); 
                Diff(i,j) = sum(sum(D));
                S(i,j)=exp(-Diff(i,j)/(2*sigm^2));
             end
        end   
         Pr(x,y,1:i,1:j)=S(1:i,1:j)/sum(S(:));  
     end
end
P=Pr;

%% compute K(x,y) = G_sigma(x-y) d(Px,Py)
fprintf('computing Motion-3d distance......');

nl_dis=zeros(nx,ny,q1,q1);
for x=1:nx
     for y=1:ny
            selx = x-q:x+q;
            sely = y-q:y+q;            

            selx(selx<1 | selx>nx | selx+dx<1 | selx+dx>nx) = [];
            sely(sely<1 | sely>ny | sely+dy<1 | sely+dy>ny) = [];
           
          %% motion distance d(Pr(x), Pr(y))
            D = 1-sum((sqrt(Pr(selx,sely,:)).*sqrt(repmat(Pr(x,y,:),[length(selx) length(sely) 1]))), 3); 
           
          %% K(x,y) = G_sigma(x-y)*d(Pr(x), Pr(y))
            K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1); 
            nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1)=K.*D; 
     end
end