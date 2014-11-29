%%%% load first frame %%%%%%
I= imread('./img/img100.png');
Img= double(I);


%%% parameter setting
timestep=5;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
iter_inner=5;
iter_outer=40;
lambda=0.5; % coefficient of the weighted length term L(phi)
strong_alpha = 0.75;  % coefficient of the weighted area term A(phi)
alfa=strong_alpha; 
epsilon=1; % papramater that specifies the width of the DiracDelta function
potentialFunction = 'double-well';
c0=2; %%%% value in the step function of initial contour 
expand_size = 10;
se_exp = strel('disk',expand_size);

%%%%%% smooth operator %%%%
sigma=1.5;     % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma);
Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.

%%%% manual initialization %%%%%
% initialize LSF as binary step function
initialLSF=c0*ones(size(Img));
initialLSF(25:100, 80:120)=-c0;  
phi=initialLSF;

%figure(2);
%imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
%title('Initial zero level contour');


% start level set evolution
for n=1:iter_outer
    phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
    %if mod(n,2)==0
    %    figure(2);
    %    imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
    %end
end


% refine the zero level contour by further level set evolution with alfa=0
alfa=0; %%%
iter_refine = 10;
phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);

%finalLSF=phi;

h=figure(1); imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); 
hold on;  contour(phi, [0,0], 'r'); hold off
%str=['Final zero level contour, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
%title(str);
saveas(h,'./seg/img100.png','png');

for i=1:1:44
    I= imread(['./img/img',num2str(100+i),'.png']);
    Img= double(I);
    Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
    [Ix,Iy]=gradient(Img_smooth);
    f=Ix.^2+Iy.^2;
    g=1./(1+f);  % edge indicator function.
    
    %%%% contour propgation %%%%
    tmp = imcomplement(phi);
    tmp = im2bw(tmp,graythresh(tmp));
    tmp = imdilate(tmp,se_exp);
    
    initialLSF=c0*ones(size(Img));
    initialLSF(tmp>0)=-c0;
    phi=initialLSF;
    
    %%%% start level set evolution %%%%
    alfa=strong_alpha;  % coefficient of the weighted area term A(phi)
    for n=1:iter_outer
        phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
    end
    
    % refine the zero level contour by further level set evolution with alfa=0
    alfa=0;  
    phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
    
    h=figure(i+1); imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); 
    hold on;  contour(phi, [0,0], 'r'); hold off
    
    saveas(h,['./seg/img',num2str(100+i),'.png'],'png');
end