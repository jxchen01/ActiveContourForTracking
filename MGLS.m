function [phi, kai] = MGLS(initial_contour, id_map, Raw, MatchingForce,betta)

%%%% parameters %%%%
sigma=1.5;     % scale parameter in Gaussian kernel
G=fspecial('gaussian',9,sigma);
c0=3;
MaxIter=500;
max_iter_inner=5;
numObj=max(id_map(:));

timestep=5;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
alfa=-1;
lambda=0.1;
epsilon=2; % papramater that specifies the width of the DiracDelta function
%betta=0;
smallNumber = 1e-10;

[dimx,dimy]=size(Raw);


%%%% edge indicator (gradient)  %%%%%
% Img_smooth=conv2(Raw,G,'same');  % smooth image by Gaussiin convolution
% [Ix,Iy]=gradient(Img_smooth);
% f=Ix.^2+Iy.^2;
% g=1./(1+f);  % edge indicator function.

opt_frangi=struct('FrangiScaleRange', [1 3], 'FrangiScaleRatio', 1, 'FrangiBetaOne',...
    0.5, 'FrangiBetaTwo', 9, 'verbose',false,'BlackWhite',false);
g=FrangiFilter2D(double(Raw),opt_frangi);
g=mat2gray(g);

%%%% intial contour %%%%%
initialLSF=c0*ones(dimx,dimy);
initialLSF(initial_contour>0)=-c0;  
phi=initialLSF;

kai = id_map;

%imagesc(Raw,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');

%%%% non-PDE version: With Matching Term %%%%
% para = struct('Raw',Raw,'iter_max',MaxIter,'mu',mu,'lambda',lambda,...
%     'alfa',alfa,'beta',betta,'epsilon',epsilon,'xdim',dimx,'ydim',dimy,...
%     'numObj',numObj,'smallNumber',smallNumber);
% [phi,kai]=LSF_fast(phi, kai, g, MatchingForce,para);

%%%% PDE Version: Only Topology-Preserving, No Matching %%%%
for iter=1:1:MaxIter
    [phi,kai] = LSF_update(phi, kai, g,lambda, mu, alfa, epsilon, timestep, max_iter_inner);
    figure(1);
    imagesc(Raw,[0, 255]); axis off; axis equal; colormap(gray); 
    hold on;  
    contour(phi, [0,0], 'r');
    drawnow
    if(iter>90 && mod(iter,10)==4)
        keyboard;
    end
end

