function nl_dis=compute_l2_dist_gray(M,nx,ny,w1,w,gau_a,q,q1,GauK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute L2 distance for gray-scaled image
%
% d(Px,Py) = int G_a(t)||Px(t) - Py(t)||^2 dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('computing L2-1d distance......');

%% extract patches
P = extract_patch(M,nx,ny,w1,w);

%% Gaussian function G_a(t) in weighted L2-distance
Gau = zeros(nx, ny, w1, w1);
Gaussi=fspecial('gaussian',[w1 w1],gau_a);
Gaussi=Gaussi/max(Gaussi(:));
for i=1:nx 
    for j=1:ny 
        Gau(i,j,1:w1,1:w1)=Gaussi(1:w1,1:w1); 
    end; 
end;

%% compute K(x,y) = G_sigma(x,y) * dist(Px,Py)
nl_dis=zeros(nx,ny,q1,q1);
for x=1:nx        
    for y=1:ny

            selx = x-q:x+q;
            sely = y-q:y+q;
            
            selx(selx<1 | selx>nx) = [];
            sely(sely<1 | sely>ny) = [];
              
          %% weighted L2-distance d(Px,Py)
            D = sum( Gau(selx,sely,:).*(P(selx,sely,:) - repmat(P(x,y,:), [length(selx) length(sely) 1])).^2, 3 );
%             D = sum( (P(selx,sely,:) - repmat(P(x,y,:), [length(selx) length(sely) 1])).^2, 3 ); %

          %% K(x,y) = G_sigma(x-y)*d(Px,Py)
            K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1);    
            nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1)=K.*D;
    end     
end

