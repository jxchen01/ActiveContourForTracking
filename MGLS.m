function id_map = MGLS(initial_contour, Raw, MatchingForce, MaxIter)

%%%% parameters %%%%
sigma=1.5;     % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma);
c0=2;

[dimx,dimy]=size(Raw);


%%%% term for gradient  %%%%%
Img_smooth=conv2(Raw,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.


%%%% intial contour %%%%%
initialLSF=c0*ones(dimx,dimy);
initialLSF(initial_contour>0)=-c0;  
phi=initialLSF;

