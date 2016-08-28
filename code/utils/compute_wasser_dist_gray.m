function nl_dis=compute_wasser_dist_gray(M,nx,ny,w1,w,q,q1,GauK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute 1d-Wasserstein (L2) distance
%
% d(Px,Py) = || sort(Px) - sort(Py) ||^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('computing wasserstein-1d distance......');

%% extract patches
P = extract_patch(M,nx,ny,w1,w);

%% Sorted vector "P"
for i=1:nx        
    for j=1:ny
              Pat(1:w1,1:w1)=P(i,j,1:w1,1:w1);
              Pat_tmp=Pat;
              Pat_sort(1:w1*w1,1) = sort(Pat_tmp(:));
              Pnew(i,j,1:w1*w1) = Pat_sort(1:w1*w1,1);        
    end
 end      
 P=Pnew;

%% compute K(x,y)=G_sigma(x,y) * dist(Px,Py)
nl_dis=zeros(nx,ny,q1,q1);
for x=1:nx        
    for y=1:ny

            selx = x-q:x+q;
            sely = y-q:y+q;
            
            selx(selx<1 | selx>nx) = [];
            sely(sely<1 | sely>ny) = [];
              
          %% d(Px,Py)
           D = sum( (P(selx,sely,:) - repmat(P(x,y,:), [length(selx) length(sely) 1])).^2, 3 );

          %% K(x,y) = G_sigma(x-y)*d(px,py)
            K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1);    
            nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1)=K.*D;
    end     
end

