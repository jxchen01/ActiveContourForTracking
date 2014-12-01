function [L_in,L_out,phi_0]=Initialization(kai,phi_0,para)

xdim=para.xdim;
ydim=para.ydim;
numObj=para.numObj;

nbr4 = [1,0; -1,0; 0,1; 0,-1];

L_in=cell(1,numObj);
for kk=1:1:numObj
    % initialize L_in for object kk
    idx_in = find(kai==kk);
    [x_in, y_in]=ind2sub([xdim,ydim], idx_in);
    head=dlnode(0);
    p=head;
    for i=1:1:numel(idx_in)
        tx=x_in(i); ty=y_in(i);
        flag=0;
        for j=1:1:4
            nx=tx+nbr4(j,1); ny=ty+nbr4(j,2);
            if(phi_0(nx,ny)>0)
                flag=1;
                break;
            end
        end
        
        if(flag)
            new_node = dlnode([tx,ty]);
            new_node.insertAfter(p);
            p=new_node;
            phi_0(tx,ty)=-1;
        else
            phi_0(tx,ty)=-3;
        end
    end
    
    tail=dlnode(0);
    tail.insertAfter(p);
    
    L_in{kk}=struct('head',head,'tail',tail);
end
    
L_out=cell(1,numObj);
for kk=1:1:numObj
    head=dlnode(0);
    tail=dlnode(0);
    tail.insertAfter(head);
    L_out{kk}=struct('head',head,'tail',tail);
end


% initialize L_out 
idx_bg=find(kai==0);
[x_bg,y_bg]=ind2sub([xdim,ydim],idx_bg);
for i=1:1:numel(idx_bg)
    tx=x_bg(i); ty=y_bg(i);
    flag=0;
    for j=1:1:4
        nx=tx+nbr4(j,1); ny=ty+nbr4(j,2);
        if(nx<1 || nx>xdim || ny<1 || ny>ydim)
            continue;
        end
        if(phi_0(nx,ny)<0)
            flag=kai(nx,ny);
            if(flag<1)
                disp('error in initializing L_out');
                keyboard;
            end
            break;
        end
    end
    
    if(flag)
        new_node = dlnode([tx,ty]);
        new_node.insertBefore(L_out{flag}.tail);
        phi_0(tx,ty)=1;
    end
end