function P = extract_patch(M,nx,ny,w1,w)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract patches P from a given (1d) image M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% location of pixels
[Y,X] = meshgrid(1:ny,1:nx);
% offsets
[dY,dX] = meshgrid(-w:w,-w:w);

% location of pixels to extract
dX = reshape(dX, [1 1 w1 w1]);
dY = reshape(dY, [1 1 w1 w1]);
X = repmat(X, [1 1 w1 w1]) + repmat(dX, [nx ny 1 1]);
Y = repmat(Y, [1 1 w1 w1]) + repmat(dY, [nx ny 1 1]);

X(X<1) = 2-X(X<1);
Y(Y<1) = 2-Y(Y<1);
X(X>ny) = 2*ny-X(X>ny);
Y(Y>nx) = 2*nx-Y(Y>nx);

% extracted patches from the image M
P = M(X + (Y-1)*ny);

