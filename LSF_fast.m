function [phi,kai]=LSF_fast(phi_0, kai, g, MF,para)

%%%% parameters %%%%
iter_max = para.iter_max;
numObj = para.numObj;
%Umax = 3;
%Vmax = 1;

%nbr4 = [1,1; 1,-1; -1,1; -1,-1];
%nbr=[1,1; 1,0; 1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];



%%%% initialization %%%%
[L_in,L_out,phi_0]=Initialization(kai,phi_0,para);

figure(2);
imagesc(para.Raw,[0, 255]); axis off; axis equal; colormap(gray); 
hold on;  
contour(phi_0, [-1,-1], 'r');
drawnow

%%%%% main loop %%%%%
phi=phi_0;
stop_flag=false;
for iter_num=1:1:iter_max
  
    Fd = computeForce(phi, g, para);
    %%%% cycle one %%%%
    
    %%%%% check L_out  %%%%%%
    [L_in,L_out,phi,kai]=check_Lout(L_in,L_out,phi,kai,Fd,MF,numObj);
    
    %%%%% update L_in %%%%%%%
    [L_in, phi]=update_Lin(L_in,phi,numObj);
    
    %%%%% check L_in %%%%%%
    [L_in,L_out,phi,kai]=check_Lin(L_in,L_out,phi,kai,Fd,MF,numObj);
    
    %%%%% update L_out %%%%%%
    [L_out, phi]=update_Lout(L_out,phi,numObj);
    
    %%%%% check stop condition %%%%
    if(stop_check(L_in,L_out,Fd,numObj))
        stop_flag=true;
    end
    
    %%%%% cycle two %%%%%%
    for repeat_num=1:1:2
        [phi,kai,L_in,L_out] = regulation(phi,kai,L_in,L_out,numObj);
    end
    
    
    figure(1);
    imagesc(para.Raw,[0, 255]); axis off; axis equal; colormap(gray); 
    hold on;  
    contour(phi, [-1,-1], 'r');
    title([num2str(iter_num),'  iteration']);
    drawnow
    
   % keyboard
    if(stop_flag)
        break;
    end
end

end


function [phi,kai,L_in,L_out] = regulation(phi,kai,L_in,L_out,numObj)
    nbr4 = [1,0; -1,0; 0,1; 0,-1];
    for obj=1:1:numObj
        p=L_out{obj}.head; OutListTail=L_out{obj}.tail;
        while(~isempty(p.Next) && length(p.Next.Data)>1)
            p=p.Next;
            nn=p.Data;
            tx=nn(1); ty=nn(2);
            if(G_phi(tx,ty,phi)<0)
                %%% switch in %%%
                if(digitTopology(tx,ty,kai,obj)>1)
                    continue;
                else                    
                    removeNode(p);
                    p.insertBefore(L_in{obj}.tail);
                    phi(tx,ty)=-1;
                    kai(tx,ty)=obj;
            
                    for i=1:1:4
                        nx=tx+nbr4(i,1); ny=ty+nbr4(i,2);
                        if(phi(nx,ny)==3)
                            phi(nx,ny)=1;
                            new_node=dlnode([nx,ny]);
                            new_node.insertBefore(OutListTail);
                        end
                    end
                end
            end
        end
    end
    
    [L_in, phi]=update_Lin(L_in,phi,numObj);
    
    for obj=1:1:numObj
        p=L_out{obj}.head; InListTail=L_out{obj}.tail;
        while(~isempty(p.Next) && length(p.Next.Data)>1)
            p=p.Next;
            nn=p.Data;
            tx=nn(1); ty=nn(2);
            if(G_phi(tx,ty,phi)>0)
                %%% switch out %%%
                if(digitTopology(tx,ty,kai,obj)>1)
                    continue;
                else                    
                    removeNode(p);
                    p.insertBefore(L_in{obj}.tail);
                    phi(tx,ty)=-1;
                    kai(tx,ty)=obj;
            
                    for i=1:1:4
                        nx=tx+nbr4(i,1); ny=ty+nbr4(i,2);
                        if(phi(nx,ny)==3)
                            phi(nx,ny)=1;
                            new_node=dlnode([nx,ny]);
                            new_node.insertBefore(InListTail);
                        end
                    end
                end
            end
        end
    end
    
    [L_out, phi]=update_Lout(L_out,phi,numObj);

end

function v= G_phi(tx,ty,phi)
    G=[1,2,1,2,4,2,1,2,1];
    nbr9=[-1,-1; -1,0; -1,1;0,-1;0,0;0,1;1,-1;1,0;1,1];
    v=0;
    for i=1:1:9
        v=v+phi(tx+nbr9(i,1), ty+nbr9(i,2))*G(i);
    end
end

function Tr=digitTopology(tx,ty,kai, objIdx)
    nbr=[1,1; 1,0; 1,-1;0,1;0,-1;-1,1;-1,0;-1,-1];
    Tn=[]; T0=0; Tobj=0;
    for j=1:1:8
        nx=tx+nbr(j,1);ny=ty+nbr(j,2);
        if(kai(nx,ny)>0)
            Tn=cat(1,Tn,kai(nx,ny));
            if(kai(nx,ny)==objIdx)
                Tobj=Tobj+1;
            end
        else
            T0=T0+1;
        end
    end
    idx_overlap=unique(Tn);
    On=numel(idx_overlap);

    Tr=min([On, max([Tobj,T0])]);
end


function [L_in,L_out,phi,kai]=check_Lin(L_in,L_out,phi,kai,Fd,MF,numObj)
    nbr4 = [1,0; -1,0; 0,1; 0,-1];
    for obj=1:1:numObj
        p=L_in{obj}.head; InListTail=L_in{obj}.tail;
        while(~isempty(p.Next) && length(p.Next.Data)>1)
            p=p.Next;
            nn=p.Data;
            tx=nn(1); ty=nn(2);
            Fm=sum(MF(tx,ty,:)) - MF(tx,ty,obj);
            if(Fd(tx,ty)-10*Fm<0)
                %%% switch out %%%
                removeNode(p);
                p.insertBefore(L_out{obj}.tail);
                phi(tx,ty)=1;
                kai(tx,ty)=0;
            
                for i=1:1:4
                    nx=tx+nbr4(i,1); ny=ty+nbr4(i,2);
                    if(phi(nx,ny)==-3)
                        phi(nx,ny)=-1;
                        new_node=dlnode([nx,ny]);
                        new_node.insertBefore(InListTail);
                    end
                end            
            end        
        end
    end
end

function [L_in,L_out,phi,kai]=check_Lout(L_in,L_out,phi,kai,Fd,MF,numObj)
    nbr4 = [1,0; -1,0; 0,1; 0,-1];
    for obj=1:1:numObj
        p=L_out{obj}.head; OutListTail=L_out{obj}.tail;
        while(~isempty(p.Next) && length(p.Next.Data)>1)
            p=p.Next;
            nn=p.Data;
            tx=nn(1); ty=nn(2);
            Fm=sum(MF(tx,ty,:)) - MF(tx,ty,obj);
            if(Fd(tx,ty)-10*Fm>0)
                %%% switch in %%%
                if(digitTopology(tx,ty,kai,obj)>1)
                    continue;
                else                    
                    removeNode(p);
                    p.insertBefore(L_in{obj}.tail);
                    phi(tx,ty)=-1;
                    kai(tx,ty)=obj;
            
                    for i=1:1:4
                        nx=tx+nbr4(i,1); ny=ty+nbr4(i,2);
                        if(phi(nx,ny)==3)
                            phi(nx,ny)=1;
                            new_node=dlnode([nx,ny]);
                            new_node.insertBefore(OutListTail);
                        end
                    end
                end
            end
        end
    end
end

function [L_in,phi]=update_Lin(L_in,phi,numObj)
    nbr4 = [1,0; -1,0; 0,1; 0,-1];
    for obj=1:1:numObj
        p=L_in{obj}.head;
        while(~isempty(p.Next) && length(p.Next.Data)>1)
            p=p.Next;
            nn=p.Data;
            tx=nn(1);ty=nn(2);
            flag=1;
            for i=1:1:4
                if(phi(tx+nbr4(i,1),ty+nbr4(i,2))>0)
                    flag=0;
                    break;
                end
            end
            if(flag)
                q=p.Next;
                removeNode(p);
                p=q;
                phi(tx,ty)=-3;
            end
        end
    end
end


function [L_out,phi]=update_Lout(L_out,phi,numObj)
    nbr4 = [1,0; -1,0; 0,1; 0,-1];
    for obj=1:1:numObj
        p=L_out{obj}.head;
        while(~isempty(p.Next) && length(p.Next.Data)>1)
            p=p.Next;
            nn=p.Data;
            tx=nn(1);ty=nn(2);
            flag=1;
            for i=1:1:4
                if(phi(tx+nbr4(i,1),ty+nbr4(i,2))<0)
                    flag=0;
                    break;
                end
            end
            if(flag)
                q=p.Next;
                removeNode(p);
                p=q;
                phi(tx,ty)=3;
            end
        end
    end
end

function flag = stop_check(L_in, L_out, Fd, numObj)
    flag=true;
    for i=1:1:numObj
        p=L_in{i}.head;
        while(~isempty(p.Next) && length(p.Next.Data)>1)
            p=p.Next;
            nn=p.Data;
            tx=nn(1); ty=nn(2);
            if(Fd(tx,ty)<0)
                flag=false;
                break;
            end
        end
        
        p=L_out{i}.head;
        while(~isempty(p.Next) && length(p.Next.Data)>1)
            p=p.Next;
            nn=p.Data;
            tx=nn(1); ty=nn(2);
            if(Fd(tx,ty)>0)
                flag=false;
                break;
            end
        end
    end
end