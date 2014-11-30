function phi = MGLS(initial_contour, Raw, MatchingForce,betta)

%%%% parameters %%%%
sigma=1.5;     % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma);
c0=2;
MaxIter=40;
max_iter_inner=5;

timestep=5;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
lambda=0.75; 
epsilon=1; % papramater that specifies the width of the DiracDelta function
%betta=0;

[dimx,dimy]=size(Raw);


%%%% edge indicator (gradient)  %%%%%
Img_smooth=conv2(Raw,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.


%%%% intial contour %%%%%
initialLSF=c0*ones(dimx,dimy);
initialLSF(initial_contour>0)=-c0;  
phi=initialLSF;

for iter=1:1:MaxIter
    phi = LSF_update(phi, g, MatchingForce,lambda, mu, alfa, betta, epsilon, timestep, max_iter_inner);
    figure(10);
    imagesc(Raw,[0, 255]); axis off; axis equal; colormap(gray); 
    hold on;  
    contour(phi, [0,0], 'r');
    hold off
end

