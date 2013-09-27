function nl_dis=compute_wasser_dist_color(M,nx,ny,w1,w,q,q1,GauK)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Sliced-Wasserstein (L2) distance
%
% d(Px,Py) = sum_{V}|| sort(Px.*V) -  sort(Py.*V) ||^2
% where V is a unit vector in R^3, and 
% thus P.*V is 1d-projection of P to the vector V
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('computing Sliced Wasserstein (3-d) distance......');

%% extract patches
P1 = extract_patch(M(:,:,1),nx,ny,w1,w);
P2 = extract_patch(M(:,:,2),nx,ny,w1,w);
P3 = extract_patch(M(:,:,3),nx,ny,w1,w);

%% sorted vector "P" (for Sliced Wasserstein distance)
% one can choose other unit vectors for projections!
%====== Random-directions for theta
%         a1=abs(randn(3,1)); a1=a1./norm(a1); 
%         a2=abs(randn(3,1)); a2=a2./norm(a2); 
%         a3=abs(randn(3,1)); a3=a3./norm(a3); 
%         a4=abs(randn(3,1)); a4=a4./norm(a4); 
%         a5=abs(randn(3,1)); a5=a5./norm(a5); 
%         a6=abs(randn(3,1)); a6=a6./norm(a6); 
%         a7=abs(randn(3,1)); a7=a7./norm(a7); 
%======== Fixed directions for theta 
        a1=[1 0 0]; a2=[0 1 0]; a3=[0 0 1];
        
for i=1:nx            
    for j=1:ny                 
        Pat1(1:w1,1:w1)=P1(i,j,1:w1,1:w1);
        Pat2(1:w1,1:w1)=P2(i,j,1:w1,1:w1);
        Pat3(1:w1,1:w1)=P3(i,j,1:w1,1:w1);
                  
        Pat_tmp1=a1(1)*Pat1+a1(2)*Pat2+a1(3)*Pat3;                
        Pat_tmp2=a2(1)*Pat1+a2(2)*Pat2+a2(3)*Pat3;                  
        Pat_tmp3=a3(1)*Pat1+a3(2)*Pat2+a3(3)*Pat3;    
%               Pat_tmp4=a4(1)*Pat1+a4(2)*Pat2+a4(3)*Pat3;
%               Pat_tmp5=a5(1)*Pat1+a5(2)*Pat2+a5(3)*Pat3;
%               Pat_tmp6=a6(1)*Pat1+a6(2)*Pat2+a6(3)*Pat3;
%               Pat_tmp7=a7(1)*Pat1+a7(2)*Pat2+a7(3)*Pat3;
                               
        Pat_sort1(1:w1*w1,1) = sort(Pat_tmp1(:));
        Pat_sort2(1:w1*w1,1) = sort(Pat_tmp2(:));
        Pat_sort3(1:w1*w1,1) = sort(Pat_tmp3(:));   
%               Pat_sort4(1:w1*w1,1) = sort(Pat_tmp4(:));
%               Pat_sort5(1:w1*w1,1) = sort(Pat_tmp5(:));
%               Pat_sort6(1:w1*w1,1) = sort(Pat_tmp6(:));
%               Pat_sort7(1:w1*w1,1) = sort(Pat_tmp7(:));
                  
        Pnew1(i,j,1:w1*w1) = Pat_sort1(1:w1*w1,1);
        Pnew2(i,j,1:w1*w1) = Pat_sort2(1:w1*w1,1);
        Pnew3(i,j,1:w1*w1) = Pat_sort3(1:w1*w1,1);                               
%               Pnew4(i,j,1:w1*w1) = Pat_sort4(1:w1*w1,1);
%               Pnew5(i,j,1:w1*w1) = Pat_sort5(1:w1*w1,1);
%               Pnew6(i,j,1:w1*w1) = Pat_sort6(1:w1*w1,1);              
%               Pnew7(i,j,1:w1*w1) = Pat_sort7(1:w1*w1,1);
            
    end    
end
P1=Pnew1;P2=Pnew2;P3=Pnew3;       
%P4=Pnew4;P5=Pnew5;P6=Pnew6;P7=Pnew7;

%% compute K(x,y) = G_sigma(x,y)*dist(Px,Py);
nl_dis=zeros(nx,ny,q1,q1);
for x=1:nx
       for y=1:ny

            selx = x-q:x+q;
            sely = y-q:y+q;
            
            selx(selx<1 | selx>nx) = [];
            sely(sely<1 | sely>ny) = [];
            
          %% d(Px,Py)
            D1 = (P1(selx,sely,:) - repmat(P1(x,y,:), [length(selx) length(sely) 1])).^2 ;              
            D2 = (P2(selx,sely,:) - repmat(P2(x,y,:), [length(selx) length(sely) 1])).^2;                 
            D3 = (P3(selx,sely,:) - repmat(P3(x,y,:), [length(selx) length(sely) 1])).^2 ;  
            D = sum( (D1+D2+D3), 3);  
%             D4 = (P4(selx,sely,:) - repmat(P4(x,y,:), [length(selx) length(sely) 1])).^2 ;                
%             D5 = (P5(selx,sely,:) - repmat(P5(x,y,:), [length(selx) length(sely) 1])).^2;                 
%             D6 = (P6(selx,sely,:) - repmat(P6(x,y,:), [length(selx) length(sely) 1])).^2 ;   
%             D7 = (P7(selx,sely,:) - repmat(P7(x,y,:), [length(selx) length(sely) 1])).^2 ;                         
%             D = sum( (D1+D2+D3+D4+D5+D6+D7), 3);
               
           %% K(x,y) = G_sigma(x,y)*dist(Px,Py);
             K = GauK(q1+1-length(selx):q1,q1+1-length(sely):q1); 
             nl_dis(x,y,q1+1-length(selx):q1,q1+1-length(sely):q1)=K.*D;       
       end       
end

