function Mg=compute_gabor_transforms(M,nx,ny,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Gabor transformed images Mg from a given image M
% 
% parameters: " sigma_start, [freq_start, freq_end] "
% "ntheta, nsigma, nfreq": # of theta (orientaion), sigma (scale), frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% fixed parameters
options.gabor_mode = 'radial';
options.gabor_mode = 'oriented';
options.add_spacial = 0;
options.iscomplex = 1;

%% choose parameters
% options.ntheta = 1;
% options.nsigma = 2;
% options.nfreq = 4;

%% 
% compute Gabor transforms (Gabor features)
fprintf('--> computing Gabor transforms: ');
[E,F] = compute_gabor_features(M,options);
N=options.ntheta*options.nsigma*options.nfreq;

%% 
% plot of Gabor features
% Ns=round(sqrt(N))+1; 
% figure;
% for i=1:N
%     subplot(Ns,Ns,i);imagesc(E(:,:,i));colormap(gray);axis image off;
% end;

%% 
% Gabor features (rescaled between 0 and 1)
Mg=zeros(nx,ny,N);
for j=1:N
    Mg(:,:,j)=rescale(E(:,:,j));
end
