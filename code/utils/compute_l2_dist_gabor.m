function nl_dis=compute_l2_dist_gabor(Mg,nx,ny,w1,w,gau_a,q,q1,GauK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute L2 distance for Gabor transformed images
%
% d(Px,Py) = int G_a(t)|| Px(t) - Py(t) ||^2 dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('computing L2-gabor distance.......');
       
%% extract Patches
P1 = extract_patch(Mg(:,:,1),nx,ny,w1,w);
P2 = extract_patch(Mg(:,:,2),nx,ny,w1,w);
P3 = extract_patch(Mg(:,:,3),nx,ny,w1,w);
P4 = extract_patch(Mg(:,:,4),nx,ny,w1,w);
% P5 = extract_patch(Mg(:,:,5),nx,ny,w1,w);
% P6 = extract_patch(Mg(:,:,6),nx,ny,w1,w);
% P7 = extract_patch(Mg(:,:,7),nx,ny,w1,w);
% P8 = extract_patch(Mg(:,:,8),nx,ny,w1,w);

%% Gaussian function G_a(t) in weighted L2-distance
Gau = zeros(nx, ny, w1, w1);
Gaussi=fspecial('gaussian',[w1 w1],gau_a);
Gaussi=Gaussi/max(Gaussi(:));
for i=1:nx 
    for j=1:ny 
        Gau(i,j,1:w1,1:w1) = Gaussi(1:w1,1:w1); 
    end
end

%% compute K(x,y) = G_sigma(x,y)*d(Px,Py)
nl_dis=zeros(nx,ny,q1,q1);
  for x=1:nx
        for y=1:ny
            selx = x-q:x+q;
            sely = y-q:y+q;
            
            selx(selx<1 | selx>nx) = [];
            sely(sely<1 | sely>ny) = [];
            
          %% weighted L2-distance d(Px,Py)
            D1 = (P1(selx,sely,:) - repmat(P1(x,y,:), [length(selx) length(sely) 1])).^2 ;                
            D2 = (P2(selx,sely,:) - repmat(P2(x,y,:), [length(selx) length(sely) 1])).^2;                 
            D3 = (P3(selx,sely,:) - repmat(P3(x,y,:), [length(selx) length(sely) 1])).^2 ;                
            D4 = (P4(selx,sely,:) - repmat(P4(x,y,:), [length(selx) length(sely) 1])).^2 ; 
            D = sum( Gau(selx,sely,:).*(D1+D2+D3+D4), 3); 
%             D5 = (P5(selx,sely,:) - repmat(P5(x,y,:), [length(selx) length(sely) 1])).^2 ;
%             D6 = (P6(selx,sely,:) - repmat(P6(x,y,:), [length(selx) length(sely) 1])).^2 ;
%             D7 = (P7(selx,sely,:) - repmat(P7(x,y,:), [length(selx) length(sely) 1])).^2 ;
%             D8 = (P8(selx,sely,:) - repmat(P8(x,y,:), [length(selx) length(sely) 1])).^2 ;                      
%             D = sum( Gau(selx,sely,:).*(D1+D2+D3+D4+D5+D6+D7+D8), 3);

          %% K(x,y) = G_sigma(x-y)*d(Px,Py)
            K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1);
            nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1)=K.*D;
        end       
  end

