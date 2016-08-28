function [M M2 phi phi2 phi3]= load_image_and_initLSF(example) 

%% choose examples 

switch example
    case 'l2-1d'         
       n=100;
       M=1+zeros(n,n);                                  
       r1=15; c=[30,30]; r2=7;  c1=[75,50]; r3=10; 
        for i=1:n
        for j=1:n
            if sqrt((i-c(1))^2+(j-c(2))^2)<r1 && sqrt((i-c(1))^2+(j-c(2))^2)> r2     
                M(i,j)=0.7;
            end
            if i>55 && i < 85 && j>25 && j <75
                if sqrt((i-c1(1))^2+(j-c1(2))^2) > r3 
                    M(i,j)=0.9-0.005*(j-25);     
                end
            end
        end        
        end
        for i=20:51
               M(i,i+35:86)=0.2;                    
        end
        M=rescale(M);

       [Y,X] = meshgrid(1:n,1:n);
       D0 = zeros(n,n)+Inf;
       c=[50,50]; r=30; %% example.1.
       D0 = -min( D0, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 )- r);
       M2 = zeros(n,n); phi2=zeros(n,n); phi3=zeros(n,n);
       
    case 'l2-gabor'
       M = load_image('ex_gabor'); 
       M = rescale((M(:,:,1)+M(:,:,2)+M(:,:,3))/3);
       [nx ny]=size(M(:,:,1)); n=nx; 

       [Y,X] = meshgrid(1:nx,1:nx);
       D0 = zeros(nx,nx)+Inf;
       r=52; c=[100,100];  %% GaborExample_last one
       D0 = -min( D0, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 ) - r );
       M2 = zeros(n,n); phi2=zeros(n,n); phi3=zeros(n,n);
       
    case 'l2-3d'
       M = load_image('ex_color'); 
       M = rescale(M);
       [nx ny]=size(M(:,:,1)); n=nx; 

       [Y,X] = meshgrid(1:nx,1:nx);
       D0 = zeros(nx,nx)+Inf; 
       r=32; c=[154,42]; 
       D0 = min( D0, sqrt( 0.31*(X-c(1)).^2 + (Y-c(2)).^2 )- r);
       r=30; c=[65,140]; 
       D0 = -min( D0, sqrt( (X-c(1)).^2 + 0.25*(Y-c(2)).^2 )- r);
       M2 = zeros(n,n); phi2=zeros(n,n); phi3=zeros(n,n);
       
   case 'wasserstein-1d'
        M = load_image('ex_wasser1d'); 
        M = rescale((M(:,:,1)+M(:,:,2)+M(:,:,3))/3);
        [nx ny]=size(M(:,:,1)); n=nx;
        
        [Y,X] = meshgrid(1:n,1:n); 
        D0 = zeros(n,n)+Inf;
        r=27; c=[33,95]; 
        D0 = min( D0, sqrt((X-c(1)).^2 + (Y-c(2)).^2 )- r);
        r=27; c=[100,40];
        D0 = -min( D0, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 )- r);
        M2 = zeros(n,n); phi2=zeros(n,n); phi3=zeros(n,n);
        
    case 'wasserstein-3d'
        M = load_image('ex_wasser3d'); 
        M = rescale(M);
        [nx ny]=size(M(:,:,1)); n=nx;
        den=0.3; choose_noise=1; %% 0 :salt & pepper, 1: random
        M = impulsenoise(M*256, den, choose_noise)/256; 

        [Y,X] = meshgrid(1:n,1:n);
        D0 = zeros(n,n)+Inf;
        r=25; c=[55,30]; 
        D0 = min( D0, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 )- r);
         r=25; c=[55,90];  
        D0 = min( D0, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 )- r);
        r=20; c=[55,140]; 
        D0 = min( D0, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 )- r);
        r=25; c=[130,40]; 
        D0 = min( D0, sqrt((X-c(1)).^2 + (Y-c(2)).^2 )- r);
        r=25; c=[128,90];
        D0 =min( D0, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 )- r);  
        r=25; c=[128,140]; 
        D0 = -min( D0, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 )- r);  
        M2 = zeros(n,n); phi2=zeros(n,n); phi3=zeros(n,n);
        
   case 'motion-1d'
        n=100;
        M=1+zeros(n,n);
        for i=1:25
            M(4*i-3:4*i-2,:) = 0;
        end
        M(30:69,30:69)=1;
        for i=1:10
            M(30:69,30+4*i-4:30+4*i-3)= 0;
        end
        u=2; v=3; 
        M2=1+zeros(n,n);
        for i=1:25
            M2(4*i-3:4*i-2,:) = 0;
        end
        M2(30+u:69+u,30+v:69+v)=1;
        for i=1:10
            M2(30+u:69+u,30+v+4*i-4:30+v+4*i-3)= 0;
        end        
        [Y,X] = meshgrid(1:n,1:n);
        D0 = zeros(n,n)+Inf;
        r=20; c=[40,40];
        D0 = -min( D0, sqrt((X-c(1)).^2 +  (Y-c(2)).^2 )- r);  
        phi2=zeros(n,n); phi3=zeros(n,n);
                           
    case 'multiphase-1d-mr'
            n=100;           
            M=zeros(n,n);
            for i=1:16; M(1:n,i)=0.2-(0.2/15)*(i-1);end
            M(1:n,17:33)=0.2;
            M(1:n,34:50)=0.4;
            M(1:n,51:67)=0.6;
            M(1:n,68:84)=0.8;
            for i=85:100; M(1:n,i)=1-(0.2/15)*(i-85); end
            M=M'; M=rescale(M);
            
            [Y,X] = meshgrid(1:n,1:n);
            D0 = zeros(n,n)+Inf;
            D02 = zeros(n,n)+Inf;
            c=[8,50]; r=7;
            D0 =min( D0, sqrt((X-c(1)).^2 +0.5*(Y-c(2)).^2 )- r);
            c=[42,50]; 
            D0 = min( D0, sqrt((X-c(1)).^2 +0.5*(Y-c(2)).^2 )- r);
             c=[77,50]; 
            D0= -min( D0, sqrt((X-c(1)).^2 +0.5*(Y-c(2)).^2 )- r);
            c=[25,50];
            D02= min( D02, sqrt( (X-c(1)).^2 +0.5* (Y-c(2)).^2 )- r);
            c=[59,50];
            D02= min( D02, sqrt( (X-c(1)).^2 + 0.5*(Y-c(2)).^2 )- r);
             c=[93,50];
            D02= min( D02, sqrt( (X-c(1)).^2 + 0.5*(Y-c(2)).^2 )- r);
            phi2 = -D02;   
            M2=zeros(n,n); phi3=zeros(n,n);
            
      case 'multiphase-1d-mi'
            n=100;
            [Y,X] = meshgrid(1:n,1:n);
            r=n*sqrt(2); c=[0,0];  % Smoothly varying intensity
            D0 = zeros(n,n)+Inf;
            D0 = - min( D0, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 )- r);
            M=rescale(D0);  
            c1=[50,40]; r1=23; c2=[50,60]; r2=23; 
            for i=1:n
                for j=1:n
                    if sqrt((i-c1(1)).^2+(j-c1(2)).^2)<r1
                        M(i,j)=0;
                    end
                    if sqrt((i-c2(1)).^2+(j-c2(2)).^2)<r2
                        M(i,j)=1;
                    end
                    if sqrt((i-c1(1)).^2+(j-c1(2)).^2)<r1 &&  sqrt((i-c2(1)).^2+(j-c2(2)).^2)<r2
                        M(i,j)=0.5;
                    end
                end
            end
            
            [Y,X] = meshgrid(1:n,1:n);
            D0 = zeros(n,n)+Inf;
            D02 = zeros(n,n)+Inf;
            r=30; c=[50,40]; 
            D0 = -min( D0, sqrt((X-c(1)).^2 +(Y-c(2)).^2 )- r);
            r=30; c=[50,60];
            D02= min( D02, sqrt( (X-c(1)).^2 + (Y-c(2)).^2 )- r);
            phi2 = -D02; 
            M2=zeros(n,n); phi3= zeros(n,n);
        
end
phi = D0;
