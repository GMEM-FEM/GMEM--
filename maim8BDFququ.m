function main
clf;
clear;
% sphere_shell4_4 
% curved_beam2_in 

filename='curved_beam2_in.I';      %input model name
element_sym=0;    %generate symetric of unsymetric element
fid=fopen(filename,'r');%打开文件
temp=fscanf(fid, '%d'); 
nel=temp(1);                 % number of elements
nnode=temp(2);                 % total number of nodes in system
nnel=8;                % number of nodes per element
ndof=3;                % number of dofs per node
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element 
temp= fgetl(fid); temp=fscanf(fid, '%f'); 
emodule= temp(1);             % elastic modulus
poisson= temp(2);             % Poisson's ratio

temp= fgetl(fid); temp= fgetl(fid); temp= fgetl(fid); temp= fgetl(fid); temp= fgetl(fid); temp= fgetl(fid); temp= fgetl(fid);  %not used
%nodes of each element
for i=1:nel
  temp= fgetl(fid);
  nodes(i,1:6) = strread(temp(20:end),' %d') ;
  temp= fgetl(fid);
  nodes(i,7:8) = strread(temp(8:end),' %d') ;
end
%coordinates of each node
temp= fgetl(fid);
for i=1:nnode
    temp= fgetl(fid); 
    gcoord(i,1:3)=strread(temp(20:end),' %f') ; 
end

ff=zeros(sdof,1);   
kk= spalloc(sdof,sdof,sdof*60);
disp=zeros(sdof,1);    eldisp=zeros(edof,1);  index=zeros(edof,1); kinmtx=zeros(6,edof);  matmtx=zeros(6,6);   4
kinmtxall=zeros(nel,2,2,2,6,24);
temp= fgetl(fid); temp= fgetl(fid); nn= strread(temp,' %d') ;disnn=fscanf(fid,'%d',nn);
%boundary
disnn1=fscanf(fid,'%d',nn*3);
disnn2=fscanf(fid,'%d',nn*3);   %not used
temp= fgetl(fid);  temp= fgetl(fid);  temp= fgetl(fid);    
conut=1;bcdof=0;bcval=0;
for i=1:nn
    inode=disnn(i);
    for j=1:3
        if(disnn1((i-1)*3+j)~=0&&(~any(bcdof==(inode-1)*3+j)))
            bcdof(conut)=(inode-1)*3+j;
            bcval(conut)=0;
            conut=conut+1;
        end
    end
end
%  force vector
nn= strread(temp,' %d') ;
disnn=fscanf(fid,'%d',nn);
disnn1=fscanf(fid,'%f',nn*3);
for i=1:nn
     inode=disnn(i);
     ff((inode-1)*3+1:(inode-1)*3+3)=disnn1((i-1)*3+1:(i-1)*3+3);
end
fclose(fid);clear  disnn1 nn conut  temp

Xe = [-1 ,1, 1 ,-1 ,-1,1 ,1 ,-1 ;
     -1,-1, 1 ,1, -1,  -1, 1 ,1;
     -1,-1 ,-1 ,-1 , 1, 1, 1 ,1];
Xe=Xe';
% sampling points & weights
point3=[-0.577350269189626,-0.577350269189626,-0.577350269189626;0.577350269189626,0.577350269189626,0.577350269189626];
weight3=[1,1,1;1,1,1];
x1=Xe(:,1);
y1=Xe(:,2);
z1=Xe(:,3);

matD= emodule/((1+poisson)*(1-2*poisson))* ...
[(1-poisson)  poisson  poisson   0   0    0; 
poisson  (1-poisson)   poisson   0   0    0;
poisson  poisson  (1-poisson)    0   0    0;
0    0    0    (1-2*poisson)/2   0    0;
0    0    0    0    (1-2*poisson)/2   0;
0    0    0    0    0   (1-2*poisson)/2];

for iel=1:nel           % loop for the total number of elements
    for i=1:nnel
        nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
        xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
        ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
        zcoord(i)=gcoord(nd(i),3);  % extract z value of the node
    end
    X(1:8,1)=xcoord;
    X(1:8,2)=ycoord;
    X(1:8,3)=zcoord;

    %axis a
    Ax0=X(1,:)+X(4,:)+X(5,:)+X(8,:)-X(2,:)-X(3,:)-X(6,:)-X(7,:);
    Ay0=X(3,:)+X(4,:)+X(7,:)+X(8,:)-X(1,:)-X(2,:)-X(5,:)-X(6,:);
    Az0=X(1,:)+X(2,:)+X(3,:)+X(4,:)-X(5,:)-X(6,:)-X(7,:)-X(8,:);
    Nx=Ax0/norm(Ax0); 
    Ny=Ay0/norm(Ay0);  
    Nz=Az0/norm(Az0);

    Xnew(:,1)=X(:,1)-mean(X(:,1));   
    Xnew(:,2)=X(:,2)-mean(X(:,2)); 
    Xnew(:,3)=X(:,3)-mean(X(:,3));

    U=U_matrix(poisson,Xnew,Nx,Ny,Nz);

    xcoord=Xnew(:,1);
    ycoord=Xnew(:,2);
    zcoord=Xnew(:,3);
    ke_unsym=zeros(24,24);
    ke_sym=zeros(24,24);

    for intx=1:2
    for inty=1:2
    for intz=1:2
        xyz0(1)=point3(intx,1);                  % sampling point in x-axis
        xyz0(2)=point3(inty,2);                  % sampling point in y-axis
        xyz0(3)=point3(intz,3);                  % sampling point in z-axis
        wtx=weight3(intx,1);               % weight in x-axis
        wty=weight3(inty,2);              % weight in y-axis
        wtz=weight3(intz,3);              % weight in z-axis

        % compute shape functions and derivatives at sampling point
        % shape functions
         shape(1)=0.125*(1-xyz0(1))*(1-xyz0(2))*(1-xyz0(3));shape(2)=0.125*(1+xyz0(1))*(1-xyz0(2))*(1-xyz0(3));
         shape(3)=0.125*(1+xyz0(1))*(1+xyz0(2))*(1-xyz0(3));shape(4)=0.125*(1-xyz0(1))*(1+xyz0(2))*(1-xyz0(3));
         shape(5)=0.125*(1-xyz0(1))*(1-xyz0(2))*(1+xyz0(3));shape(6)=0.125*(1+xyz0(1))*(1-xyz0(2))*(1+xyz0(3));
         shape(7)=0.125*(1+xyz0(1))*(1+xyz0(2))*(1+xyz0(3));shape(8)=0.125*(1-xyz0(1))*(1+xyz0(2))*(1+xyz0(3));

        % derivatives
         dhdr(1)=-0.125*(1-xyz0(2))*(1-xyz0(3));dhdr(2)=0.125*(1-xyz0(2))*(1-xyz0(3));
         dhdr(3)=0.125*(1+xyz0(2))*(1-xyz0(3));dhdr(4)=-0.125*(1+xyz0(2))*(1-xyz0(3));
         dhdr(5)=-0.125*(1-xyz0(2))*(1+xyz0(3));dhdr(6)=0.125*(1-xyz0(2))*(1+xyz0(3));
         dhdr(7)=0.125*(1+xyz0(2))*(1+xyz0(3));dhdr(8)=-0.125*(1+xyz0(2))*(1+xyz0(3));
         dhds(1)=-0.125*(1-xyz0(1))*(1-xyz0(3));dhds(2)=-0.125*(1+xyz0(1))*(1-xyz0(3));
         dhds(3)=0.125*(1+xyz0(1))*(1-xyz0(3));dhds(4)=0.125*(1-xyz0(1))*(1-xyz0(3));
         dhds(5)=-0.125*(1-xyz0(1))*(1+xyz0(3));dhds(6)=-0.125*(1+xyz0(1))*(1+xyz0(3));
         dhds(7)=0.125*(1+xyz0(1))*(1+xyz0(3));dhds(8)=0.125*(1-xyz0(1))*(1+xyz0(3));
         dhdt(1)=-0.125*(1-xyz0(1))*(1-xyz0(2));dhdt(2)=-0.125*(1+xyz0(1))*(1-xyz0(2));
         dhdt(3)=-0.125*(1+xyz0(1))*(1+xyz0(2));dhdt(4)=-0.125*(1-xyz0(1))*(1+xyz0(2));
         dhdt(5)=0.125*(1-xyz0(1))*(1-xyz0(2));dhdt(6)=0.125*(1+xyz0(1))*(1-xyz0(2));
         dhdt(7)=0.125*(1+xyz0(1))*(1+xyz0(2));dhdt(8)=0.125*(1-xyz0(1))*(1+xyz0(2));

        % jacob3=fejacob3(nnel,dhdr,dhds,dhdt,xcoord,ycoord,zcoord);  % compute Jacobian
         jacob3=zeros(3,3);
         for i=1:nnel
             jacob3(1,1)=jacob3(1,1)+dhdr(i)*xcoord(i);jacob3(1,2)=jacob3(1,2)+dhdr(i)*ycoord(i);jacob3(1,3)=jacob3(1,3)+dhdr(i)*zcoord(i);
             jacob3(2,1)=jacob3(2,1)+dhds(i)*xcoord(i);jacob3(2,2)=jacob3(2,2)+dhds(i)*ycoord(i);jacob3(2,3)=jacob3(2,3)+dhds(i)*zcoord(i);
             jacob3(3,1)=jacob3(3,1)+dhdt(i)*xcoord(i);jacob3(3,2)=jacob3(3,2)+dhdt(i)*ycoord(i);jacob3(3,3)=jacob3(3,3)+dhdt(i)*zcoord(i);
         end
        detjacob=det(jacob3);                 % determinant of Jacobian
        invjacob=inv(jacob3);                 % inverse of Jacobian matrix 

        for i=1:nnel
            dhdx(i)=invjacob(1,1)*dhdr(i)+invjacob(1,2)*dhds(i)+invjacob(1,3)*dhdt(i);
            dhdy(i)=invjacob(2,1)*dhdr(i)+invjacob(2,2)*dhds(i)+invjacob(2,3)*dhdt(i);
            dhdz(i)=invjacob(3,1)*dhdr(i)+invjacob(3,2)*dhds(i)+invjacob(3,3)*dhdt(i);
        end
        kinmtx1=zeros(6,24);
        % determine the kinematic equation between strains and displacements for three-dimensional solids 
         for i=1:nnel
             i1=(i-1)*3+1;  
             i2=i1+1;
             i3=i2+1;
             kinmtx1(1,i1)=dhdx(i);kinmtx1(2,i2)=dhdy(i);kinmtx1(3,i3)=dhdz(i);
             kinmtx1(4,i1)=dhdy(i);kinmtx1(4,i2)=dhdx(i);
             kinmtx1(5,i2)=dhdz(i);kinmtx1(5,i3)=dhdy(i);
             kinmtx1(6,i1)=dhdz(i);kinmtx1(6,i3)=dhdx(i);
         end

        kinmtx1=kinmtx1';
        xyz=zeros(3,1);
        for ii=1:8
            xyz(1)=xyz(1)+1/8*(1+xyz0(1)*x1(ii))*(1+xyz0(2)*y1(ii))*(1+xyz0(3)*z1(ii))* Xnew(ii,1);
            xyz(2)=xyz(2)+1/8*(1+xyz0(1)*x1(ii))*(1+xyz0(2)*y1(ii))*(1+xyz0(3)*z1(ii))* Xnew(ii,2);
            xyz(3)=xyz(3)+1/8*(1+xyz0(1)*x1(ii))*(1+xyz0(2)*y1(ii))*(1+xyz0(3)*z1(ii))* Xnew(ii,3);
        end

        for j4=3:8
            for i4=1:3   
                TTG=T_stress(xyz,poisson,i4,j4,Nx,Ny,Nz);   
                TT2(1,1)=TTG(1,1);   TT2(2,1)=TTG(2,2);   TT2(3,1)=TTG(3,3);
                TT2(4,1)=2*TTG(1,2); TT2(5,1)=2*TTG(2,3); TT2(6,1)=2*TTG(1,3);            
                kinmtx(:,(j4-1)*3+i4)=TT2(:,1);
            end
        end

        if(element_sym==1)
            ke_sym=ke_sym+wtx*wty*wtz*detjacob*kinmtx'*matD*kinmtx;    % element matrix
        else
            ke_unsym=ke_unsym+wtx*wty*wtz*detjacob*kinmtx1*matD*kinmtx;    % element matrix
        end
    end
    end
    end 
    invU=inv(U); 
    if(element_sym==1)
        k=invU'*ke_sym*invU;
    else
        k=ke_unsym*invU;
    end
    % extract system dofs associated with element
    h=0;
    for i=1:nnel
     start = (nd(i)-1)*ndof;
       for j=1:ndof
         h=h+1;
         index(h)=start+j;
       end
    end
    kk(index,index)=kk(index,index)+k;
 end
%   apply boundary conditions
 n=length(bcdof);
 sdof=size(kk);
 for i=1:n
    c=bcdof(i);
    kk(c,1:sdof)=0;
    kk(1:sdof,c)=0;
    kk(c,c)=1e8;
    ff(c)=bcval(i);
 end

%  solve the matrix equation
disp=kk\ff;
disp1=max(abs(disp(1:3:sdof))); disp2=max(abs(disp(2:3:sdof))); disp3=max(abs(disp(3:3:sdof))); 
end
%modal displacement matrix U
function [U]=U_matrix(miu,X,Nx,Ny,Nz)
U=zeros(24,24);
for k=1:3
    l=mod(k,3)+1;
    m=mod(l,3)+1; 
	for i=1:8
        %traslate
        U((i-1)*3+k,k)=1;
        %rotate
        U((i-1)*3+l,3+k)=-X(i,m);
        U((i-1)*3+m,3+k)=X(i,l);
       %tension
        U((i-1)*3+k,6+k)=X(i,k);
        U((i-1)*3+l,6+k)=-miu*X(i,l);
        U((i-1)*3+m,6+k)=-miu*X(i,m);
        %shear
        U((i-1)*3+k,15+k)=X(i,l);
    end
    % bending
    N(1,:)=Nx;N(2,:)=Ny;N(3,:)=Nz;
    T(1,:)=N(k,:);
    T(2,:)=cross(N(k,:),N(l,:));T(2,:)=T(2,:)/norm(T(2,:));
    T(3,:)=cross(N(k,:),T(2,:));
    for i=1:8
        TT((i-1)*3+1:(i-1)*3+3,(i-1)*3+1:(i-1)*3+3)=T;
        Xe_bending(i,:)=T*X(i,:)';
        U((i-1)*3+1,9+k)=-Xe_bending(i,1)*Xe_bending(i,2);
        U((i-1)*3+2,9+k)=(Xe_bending(i,1)*Xe_bending(i,1)-miu*(Xe_bending(i,3)^2-Xe_bending(i,2)^2))/2;
        U((i-1)*3+3,9+k)=miu*Xe_bending(i,3)*Xe_bending(i,2);
    end
    U(:,9+k)=TT'*U(:,9+k);
       
    T(1,:)=N(k,:);
    T(2,:)=cross(N(k,:),N(m,:));T(2,:)=T(2,:)/norm(T(2,:));
    T(3,:)=cross(N(k,:),T(2,:));
    for i=1:8
      TT((i-1)*3+1:(i-1)*3+3,(i-1)*3+1:(i-1)*3+3)=T;
      Xe_bending(i,:)=T*X(i,:)';
        U((i-1)*3+1,12+k)=-Xe_bending(i,1)*Xe_bending(i,2);
        U((i-1)*3+2,12+k)=(Xe_bending(i,1)*Xe_bending(i,1)-miu*(Xe_bending(i,3)^2-Xe_bending(i,2)^2))/2;
        U((i-1)*3+3,12+k)=miu*Xe_bending(i,3)*Xe_bending(i,2);
    end
    U(:,12+k)=TT'*U(:,12+k);
    
    % non-physical
    Ny0=cross(N(k,:),N(l,:)); Ny0=cross(Ny0,N(k,:)); Ny0=Ny0/norm(Ny0); 
    Nz0=cross(N(k,:),N(m,:));Nz0=cross(Nz0,N(k,:)); Nz0=Nz0/norm(Nz0);
    Nyz=Ny0+Nz0; Nyz=Nyz/norm(Nyz); 
    N2=cross(N(k,:),-Nyz);  N2=N2+Nyz;
    N2=N2/norm(N2); 

    N3=cross(N(k,:),Nyz);  N3=N3+Nyz;
    N3=N3/norm(N3);
    T(1,:)=N(k,:);
    T(2,:)=N2;
    T(3,:)=N3;

    for i=1:8
      TT((i-1)*3+1:(i-1)*3+3,(i-1)*3+1:(i-1)*3+3)=T;
    end
    for i=1:8
        Xe_non(i,:)=T*X(i,:)';
        U((i-1)*3+1,21+k)=Xe_non(i,1)*Xe_non(i,2)*Xe_non(i,3);
    end 

    U(:,21+k)=TT'*U(:,21+k);
    % torsion
    Xe_torsion=Xe_non;
    if(k==1||k==2)   %19~21
        for i=1:8
            U((i-1)*3+2,18+k)=Xe_torsion(i,1)*Xe_torsion(i,3);
            U((i-1)*3+3,18+k)=-Xe_torsion(i,1)*Xe_torsion(i,2);
        end
        U(:,18+k)=TT'*U(:,18+k);
    else
        for i=1:8
            U((i-1)*3+1,21)=Xe_torsion(i,2)*Xe_torsion(i,3);
            U((i-1)*3+2,21)=Xe_torsion(i,3)*Xe_torsion(i,1);
            U((i-1)*3+3,21)=Xe_torsion(i,1)*Xe_torsion(i,2);
        end
        U(:,21)=TT'*U(:,21);
    end
end            
end

function[T]=T_stress(X,miu,ii,jj,Nx,Ny,Nz)
% ii which direction, varies from 1 to 3
% jj which mode, varies from 1 to 8
l=mod(ii,3)+1;
m=mod(l,3)+1;
 N(1,:)=Nx;N(2,:)=Ny;N(3,:)=Nz;  
T=zeros(3,3);
if(jj==3)
    T(ii,ii)=1;T(l,l)=-miu;T(m,m)=-miu;
elseif(jj==4)   %xy平面
    l=mod(ii,3)+1;
    m=mod(l,3)+1;
    TT(1,:)=N(ii,:);
    TT(2,:)=cross(N(ii,:),N(l,:));TT(2,:)=TT(2,:)/norm(TT(2,:));
    TT(3,:)=cross(N(ii,:),TT(2,:));
    Xe=TT*X; 
    T(1,1)=-Xe(2); T(2,2)=miu*Xe(2);T(3,3)=miu*Xe(2);
    T=TT'*T*TT;
elseif(jj==5)  %xz平面
    l=mod(ii,3)+1;
    m=mod(l,3)+1;
    TT(1,:)=N(ii,:);
    TT(2,:)=cross(N(ii,:),N(m,:));TT(2,:)=TT(2,:)/norm(TT(2,:));
    TT(3,:)=cross(N(ii,:),TT(2,:));
    Xe=TT*X; 
    T(1,1)=-Xe(2);  T(2,2)=miu*Xe(2);T(3,3)=miu*Xe(2);
    T=TT'*T*TT;
elseif(jj==6)
    T(ii,l)=0.5;T(l,ii)=T(ii,l);
elseif(jj==7)
    l=mod(ii,3)+1;
    m=mod(l,3)+1;
    Ny0=cross(N(ii,:),N(l,:));Ny0=cross(Ny0,N(ii,:)); Ny0=Ny0/norm(Ny0); 
    Nz0=cross(N(ii,:),N(m,:));Nz0=cross(Nz0,N(ii,:));Nz0=Nz0/norm(Ny0); 
    Nyz=Ny0+Nz0; Nyz=Nyz/norm(Nyz); 
    N2=cross(N(ii,:),-Nyz);  N2=N2+Nyz;
    N2=N2/norm(N2); 
    N3=cross(N(ii,:),Nyz);  N3=N3+Nyz;
    N3=N3/norm(N3);
    TT(1,:)=N(ii,:);
    TT(2,:)=N2;
    TT(3,:)=N3;
    Xe=TT*X;
    if(ii==1)
        T(1,2)=0.5*Xe(3);T(2,1)=T(1,2);
        T(1,3)=-0.5*Xe(2);T(3,1)=T(1,3);
    elseif(ii==2)
        T(1,2)=0.5*Xe(3);T(2,1)=T(1,2);
        T(1,3)=-0.5*Xe(2);T(3,1)=T(1,3);
    else
        x=Xe(1);y=Xe(2);z=Xe(3);
        T(1,2)= z;  T(2,1)=T(1,2);
        T(2,3)= x;    T(3,2)=T(2,3);
        T(1,3)=  y;   T(3,1)=T(1,3);
    end
    T=TT'*T*TT;   
elseif(jj==8)
    l=mod(ii,3)+1;
    m=mod(l,3)+1;
    Ny0=cross(N(ii,:),N(l,:));Ny0=cross(Ny0,N(ii,:)); Ny0=Ny0/norm(Ny0); 
    Nz0=cross(N(ii,:),N(m,:));Nz0=cross(Nz0,N(ii,:)); Nz0=Nz0/norm(Nz0); 
    Nyz=Ny0+Nz0; Nyz=Nyz/norm(Nyz); 
    N2=cross(N(ii,:),-Nyz);  N2=N2+Nyz;
    N2=N2/norm(N2); 
    N3=cross(N(ii,:),Nyz);  N3=N3+Nyz;
    N3=N3/norm(N3);
    TT(1,:)=N(ii,:);
    TT(2,:)=N2;
    TT(3,:)=N3;
    Xe=TT*X;
    x=Xe(1);y=Xe(2);z=Xe(3);          
    T(1,1)= y*z;
    T(2,2)= -miu*y*z;
    T(3,3)= -miu*y*z;
    T=TT'*T*TT;  
end
end
