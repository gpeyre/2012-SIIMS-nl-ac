function [E,F] = compute_gabor_features(M,options)

% compute_gabor_features - compute a set of gabor filtered features
%
%   [E,F] = compute_gabor_features(M,options);
%
%   Each F(:,:,k) is a filtered image.
%   Each E(:,:,k) is a fature (non-linearly remapped filtering).
%
%   Copyright (c) 2007 Gabriel Peyre
%% Modified (Dec.2010) E = sqrt(abs(F).^2); 

gabor_mode =    getoptions(options, 'gabor_mode', 'oriented');
iscomplex =     getoptions(options, 'iscomplex', 1);
ntheta =        getoptions(options, 'ntheta', 8);
nsigma =        getoptions(options, 'nsigma', 6);
nfreq =         getoptions(options, 'nfreq', 4);
scaling =       getoptions(options, 'scaling',sqrt(2));
add_spacial =   getoptions(options, 'add_spacial',0);

% added:::
sigma_start = getoptions(options, 'sigma_start', '4');
freq_start = getoptions(options, 'freq_start', '2');
freq_end = getoptions(options, 'freq_end', '6');

if not(strcmp(gabor_mode,'oriented'))
    ntheta = 1;
end

n = size(M,1);
% size of the filters
m = 51; 

theta_list = linspace(0,pi,ntheta+1); theta_list(end)=[];
theta_list

% sigma_list = 2* scaling.^(0:nsigma-1); %% original
sigma_list = sigma_start* scaling.^(0:nsigma-1); 
sigma_list

if nfreq>1
%     freq_list = linspace(0.7,2.5,nfreq+1); freq_list(end)=[]; %% original
    freq_list = linspace(freq_start,freq_end,nfreq+1); freq_list(end)=[];
    freq_list
else
    freq_list = 1;
end
q = ntheta*nsigma*nfreq;
F = zeros(n,n,ntheta,nsigma,nfreq);
num = 0;
for ifreq=1:nfreq
for itheta=1:ntheta
    for isigma=1:nsigma
        num = num+1;
        progressbar(num,q);
        sigma = sigma_list(isigma);
        if strcmp(gabor_mode,'oriented')
            theta = theta_list(itheta);
        else
            theta = [];
        end
        h = compute_gabor_filter(m,sigma,theta,freq_list(ifreq)/sigma,options);
        F(:,:,itheta,isigma,ifreq) = perform_convolution(M,h,options);
    end
end
end
F = reshape(F,[n,n,q]);

add_spacial = 0;

%% compute the feature vector
% E = abs(F).^2; %% original
E = sqrt(abs(F).^2); 

if add_spacial
    [Y,X] = meshgrid(1:n,1:n);
    E(end+1:end+2,:) = [X(:),Y(:)]';
end