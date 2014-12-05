I=imread('img1_bw.png');
[dimx,dimy]=size(I);

bw0= (I>1e-8);
ctl0=bwmorph(bw0,'thin',Inf);
cc= bwconncomp(ctl0);
labmat = labelmatrix(cc);
id_map = zeros(dimx,dimy);
se=strel('disk',2);
for i=1:1:3
    im = ismember(labmat,i);
    pList = SortCellPixel(im,1);
    len = length(pList);
    sList = pList(floor(2*len/5):1:ceil(3*len/5));
    tmp=zeros(dimx,dimy);
    tmp(sList)=1;
    tmp=imdilate(tmp,se);
    id_map = id_map + tmp.*i;
end
initial_contour = (id_map>0);


Raw=imread('img2.png');

MF=zeros(dimx,dimy);
[phi,kai]=MGLS(initial_contour, id_map, Raw, MF, 0);


% MF=zeros(dimx,dimy,3);
% 
% I1=imread('img2_bw1.png');
% bw1=(I1>0);
% ctl1=bwmorph(bw1,'thin',Inf);
% dist1=bwdist(ctl1);
% dist1(dist1<1e-8)=0.01;
% tmp=1./(1+exp(-1.*(5-dist1(:,:))));
% tmp(dist1>10)=0;
% MF(:,:,1)=tmp(:,:);
% 
% 
% I2=imread('img2_bw2.png');
% bw2=(I2>0);
% ctl2=bwmorph(bw2,'thin',Inf);
% dist2=bwdist(ctl2);
% dist2(dist2<1e-8)=0.01;
% tmp=1./(1+exp(-1.*(5-dist2(:,:))));
% tmp(dist2>10)=0;
% MF(:,:,2)=tmp(:,:);
% 
% I3=imread('img2_bw3.png');
% bw3=(I3>0);
% ctl3=bwmorph(bw3,'thin',Inf);
% dist3=bwdist(ctl3);
% dist3(dist3<1e-8)=0.01;
% tmp=1./(1+exp(-1.*(5-dist3(:,:))));
% tmp(dist3>10)=0;
% MF(:,:,3)=tmp(:,:);
% 
% [phi,kai]=MGLS(initial_contour, id_map, Raw, MF, 0);



