function pixelList = SortCellPixel(im, outputType)

%%%% 
% outputType:  1 -- index,  2 -- subscription
%%%%

im(im>0)=1;
[xdim,ydim]=size(im);
nbr=[1,1;1,0;1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];

epImg = bwmorph(im, 'endpoint');
labelImg=epImg+im;

if(nnz(epImg)~=2)
    disp('error in number of endpoints');
    keyboard
end
clear epImg

ep=find(labelImg==2);


% extract a new end point
[xt,yt]=ind2sub([xdim ydim],ep(1));
labelImg(xt,yt)=0;

% extract the cell body
tmpPixInd=zeros(100,2);
pixNum=1;
tmpPixInd(1,1)=xt; tmpPixInd(1,2)=yt;
flag=1;
while(flag)
    flag=0;
    for i=1:8
        x0=xt+nbr(i,1);
        y0=yt+nbr(i,2);
        if(x0<1 || x0>xdim || y0<1 || y0>ydim)
            continue;
        end
        if(labelImg(x0,y0)>0)
            if(labelImg(x0,y0)==2)
                flag=2;
            else
                flag=1;
            end
            
            pixNum=pixNum+1;
            tmpPixInd(pixNum,1)=x0;tmpPixInd(pixNum,2)=y0;
            labelImg(x0,y0)=0;
            xt=x0; yt=y0;
            break;
        end
    end
    if(flag==2)
        break;
    end
end

if(outputType==1)
    pixelList=sub2ind([xdim,ydim],tmpPixInd(1:pixNum,1),tmpPixInd(1:pixNum,2));
elseif(outputType==2)
    pixelList=zeros(pixNum,2);
    pixelList(:,:)=tmpPixInd(1:pixNum,:);
else
    error('error output type');
end