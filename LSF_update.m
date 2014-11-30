function [phi,kai_new] = LSF_update(phi_0, kai, g, MF, lambda,mu, alfa, beta, epsilon, timestep, iter)

NB_width=timestep+4;
[dimx,dimy]=size(g);
epsilon_merge = 0.1;
nbr=[1,1; 1,0; 1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];
phi=phi_0;

[vx, vy]=gradient(g);
for k=1:iter
    phi=NeumannBoundCond(phi);
    
    [phi_x,phi_y]=gradient(phi);
    s=sqrt(phi_x.^2 + phi_y.^2);
    smallNumber=1e-10;  
    Nx=phi_x./(s+smallNumber); % add a small positive number to avoid division by zero
    Ny=phi_y./(s+smallNumber);
    curvature=div(Nx,Ny);

    distRegTerm=distReg_p2(phi);  % compute the distance regularization term in eqaution (13) with the double-well potential p2.
      
    diracPhi=Dirac(phi,epsilon);
    areaTerm=diracPhi.*g; % balloon/pressure force
    edgeTerm=diracPhi.*(vx.*Nx+vy.*Ny) + diracPhi.*g.*curvature;
    matchingTerm = diracPhi.*MF;
    phi_temp= phi + timestep*(mu*distRegTerm + lambda*edgeTerm + alfa*areaTerm + beta*matchingTerm);
    
        
    %%%% extract narrow band points %%%%
    nb = find(phi<NB_width & phi>-NB_width);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    num_nb = numel(nb);
    [nbx,nby]=ind2sub([dimx,dimy],nb);
    for i=1:1:num_nb
        tx=nbx(i);ty=nby(i);
        if(phi(tx,ty)>0 && phi_temp(tx,ty)<0)  
            Tn=[]; T0=0; 
            for j=1:1:8
                nx=tx+nbr(j,1);ny=ty+nbr(j,2);
                if(kai(nx,ny)>0)
                    Tn=cat(1,Tn,kai(nx,ny));
                else
                    T0=T0+1;
                end
            end
            idx_overlap=unique(Tn);
            On=numel(idx_overlap);                  
            if(On>1) %%% relaxed simple point
                phi(tx,ty)=epsilon_merge;
            %elseif %%% the other case is nearly impossible
            else
                phi(tx,ty)=phi_temp(tx,ty);
            end
        elseif(phi(tx,ty)<0 && phi_temp(tx,ty)>0)    
            phi(tx,ty)=phi_temp(tx,ty);
        else
            phi(tx,ty)=phi_temp(tx,ty);
        end
    end
    
    bw = (phi<1e-8);
    cc = bwconncomp(bw,4);
    kai_new=zeros(dimx,dimy);
    for i=1:1:cc.NumObjects
        idx=cc.PixelIdxList{i};
        tid=kai(idx);
        overlapID=unique(nonzeros(tid));
        numID=numel(overlapID);
        if(numID==1)
            kai_new(idx)=overlapID;
        elseif(numID>1)
            maxNum=0;
            for j=1:1:numID
                tmp=nnz(tid==overlapID(j));
                if(tmp>maxNum)
                    maxNum=tmp;
                    maxID=overlapID(j);
                end
            end
            kai_new(idx)=maxID;
        else
            disp('error in updating Kai');
            keyboard;
        end
    end
    
    kai=kai_new;
    
end


function f = distReg_p2(phi)
% compute the distance regularization term with the double-well potential p2 in eqaution (16)
[phi_x,phi_y]=gradient(phi);
s=sqrt(phi_x.^2 + phi_y.^2);
a=(s>=0) & (s<=1);
b=(s>1);
ps=a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  % compute first order derivative of the double-well potential p2 in eqaution (16)
dps=((ps~=0).*ps+(ps==0))./((s~=0).*s+(s==0));  % compute d_p(s)=p'(s)/s in equation (10). As s-->0, we have d_p(s)-->1 according to equation (18)
f = div(dps.*phi_x - phi_x, dps.*phi_y - phi_y) + 4*del2(phi);  

function f = div(nx,ny)
[nxx,~]=gradient(nx);  
[~,nyy]=gradient(ny);
f=nxx+nyy;

function f = Dirac(x, sigma)
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;

function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  

