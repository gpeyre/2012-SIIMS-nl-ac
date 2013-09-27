function Edge = edge_detector(M, compute_edge, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Edge detector (for GAC, IAC model)
% compute_edge = 'edge-1d' (gray-valued), 'edge-3d' (color), 
%                'edge-gabor' (Gabor transforms)
%
%% NOTE: gau_sigma/4 = usual Gaussian sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('computing Edge detector......');

alpha = options.alpha;
epsilon = options.epsilon;
gau_sig = options.gau_sig;

options.order=2;
switch compute_edge
    case 'edge-1d'
        G = grad(M,options);
        d = perform_blurring( sqrt(sum(G .^2,3)),gau_sig); %% case.1.
        % d = perform_blurring( (sum(G.^2,3)),gau_sig);
        Edge = (epsilon+d).^(-alpha);
        
    case 'edge-3d'
        G1 = grad(M(:,:,1),options);
        G2 = grad(M(:,:,2),options);
        G3 = grad(M(:,:,3),options);
        d = perform_blurring( sqrt(sum(G1.^2+G2.^2+G3.^2,3)),gau_sig); %% case.1.
        Edge = (epsilon+d).^(-alpha);      
        
    case 'edge-gabor'           
        [nx,ny,L]=size(M);
        d=zeros(nx,ny);
        for i=1:L
            Gradfunc= grad(perform_blurring(M(1:nx,1:ny,i),gau_sig),options); %% case.2. (this one is better)
            d = d+sum(Gradfunc.^2,3);
        end
        Edge = (epsilon+d).^(-alpha);
end
Edge = rescale(Edge); %% edge is rescaled between [0,1]
