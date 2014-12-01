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

% I1=imread('img2_bw1.png');
% bw1=(I1>0);
% ctl1=bwmorph(bw1,'thin',Inf);
% 
% I2=imread('img2_bw2.png');
% bw2=(I2>0);
% ctl2=bwmorph(bw2,'thin',Inf);






