function A=curlperiodicv2new(Nxt,Nyt,Nz,dxt,dyt,dz,kxt,kyt,kz,Mx,My,M3)

%{
close all; 
clear all; 
clc; 

Nxt=5;
Nyt=5;
Nz=1; 
dxt=0.1; 
dyt=0.1;
dz=0.1; 
kxt=0; 
kyt=0; 
kz=0; 

%we need to change definition of this if we change the way our mirror plane
%is constructed
Mx=0; 
My=1; 
M3=0; 
%}


%order of elements screwed up by other guy in other section of code, we
%should fix later

M1=My; 
M2=Mx; 
Nx=Nyt;
Ny=Nxt;
dx=dyt;
dy=dxt;
kx=kyt;
ky=kxt;

if (Nx==1)
    Dx=0;
    DxH=0; 
else
    if(M1==0)
        ox=ones(Nx,1);
        Dx = spdiags([-ox ox]/dx, 0:1,Nx,Nx);
        Dx(end,1)=1/dx*exp(1i*kx*Nx*dx);
        DxH=Dx';
    else
        ox=ones(Nx,1);
        Dx = spdiags([-ox, ox]/dx, 0:1,Nx,Nx); %differention for E;
        Dx(end,end)=2*Dx(end,end); 
        DxH = -spdiags([-ox,ox]/dx, -1:0,Nx,Nx);
    end
end

if (Ny==1)
    Dy=0;
    DyH=0; 
else
    if(M2==0)
        oy=ones(Ny,1);
        Dy = spdiags([-oy oy]/dy, 0:1,Ny,Ny);
        Dy(end,1)=1/dy*exp(1i*ky*Ny*dy);
        DyH=Dy';
    else
        oy=ones(Ny,1);
        Dy = spdiags([-oy, oy]/dy, 0:1,Ny,Ny); %differention for E;
        Dy(end,end)=2*Dy(end,end); 
        DyH = -spdiags([-oy,oy]/dy, -1:0,Ny,Ny);
    end
end

if (Nz==1)
    Dz=0;
    DzH=0; 
else
    if(M3==0)
        oz=ones(Nz,1);
        Dz = spdiags([-oz oz]/dz, 0:1,Nz,Nz);
        Dz(end,1)=1/dz*exp(1i*kz*Nz*dz);
        DzH=Dz';
    else
        oz=ones(Nz,1);
        Dz = spdiags([-oz, oz]/dz, 0:1,Nz,Nz); %differention for E;
        Dz(end,end)=2*Dz(end,end); 
        DzH = -spdiags([-oz,oz]/dz, -1:0,Nz,Nz);
    end
end

Tx = kron(kron(Dx,speye(Ny)),speye(Nz));
TxH = kron(kron(DxH,speye(Ny)),speye(Nz));

Ty = kron(kron(speye(Nx),Dy),speye(Nz));
TyH = kron(kron(speye(Nx),DyH),speye(Nz));
  
Tz = kron(speye(Nx),kron(speye(Ny),Dz));
TzH = kron(speye(Nx),kron(speye(Ny),DzH));

Nxyz= Nx*Ny*Nz;
Mzero = 0*speye(Nxyz);

CE = [Mzero,-Tz, Ty;
     Tz, Mzero, -Tx;
     -Ty,Tx,Mzero];

CH = [Mzero,TzH,-TyH;
      -TzH,Mzero, TxH;
      TyH,-TxH,Mzero];
 
% CH = [Mzero,Tz',-Ty';
%       -Tz',Mzero, Tx';
%       Ty',-Tx',Mzero];
 
%CH = CE';

Mzero2 = 0*speye(3*Nxyz);


A = [Mzero2,1i*(CH); -1i*(CE),Mzero2];

Ap=A; 

%we want to set certain components of A to zero so we can enforce the E=0,
%H=0 condition, for different boundary condition cases


if(M1==1)
    for ny=1:1:Ny
        for nz=1:1:Nz
            A(:,nz+(ny-1)*Nz+1*Nx*Ny*Nz)=0; %Ey component
            A(:,nz+(ny-1)*Nz+2*Nx*Ny*Nz)=0; %Ez component
            A(:,nz+(ny-1)*Nz+3*Nx*Ny*Nz)=0; %Hx component
        
            A(nz+(ny-1)*Nz+1*Nx*Ny*Nz,:)=0; %Ey component
            A(nz+(ny-1)*Nz+2*Nx*Ny*Nz,:)=0; %Ez component
            A(nz+(ny-1)*Nz+3*Nx*Ny*Nz,:)=0; %Hx component
        end
    end
end

if(M2==1)
    for nx=1:1:Nx
        for nz=1:1:Nz
            A(:,nz+(nx-1)*Nz*Ny)=0; %this will take care of all Ex components
            A(:,nz+(nx-1)*Nz*Ny+2*Nx*Ny*Nz)=0; %Ez
            A(:,nz+(nx-1)*Nz*Ny+4*Nx*Ny*Nz)=0;  %Hy
        
            A(nz+(nx-1)*Nz*Ny,:)=0; %this will take care of all Ex components
            A(nz+(nx-1)*Nz*Ny+2*Nx*Ny*Nz,:)=0; %Ez
            A(nz+(nx-1)*Nz*Ny+4*Nx*Ny*Nz,:)=0;  %Hy
        end
    end
end


if(M3==1)
    for ny=1:1:Ny
        for nx=1:1:Nx
            A(:,(ny-1)*Nz+(nx-1)*Nz*Ny+1)=0;  %Ex
            A(:,(ny-1)*Nz+(nx-1)*Nz*Ny+1+Nx*Ny*Nz)=0;  %Ey
            A(:,(ny-1)*Nz+(nx-1)*Nz*Ny+1+5*Nx*Ny*Nz)=0;  %Hz
        
            A((ny-1)*Nz+(nx-1)*Nz*Ny+1,:)=0;  %Ex
            A((ny-1)*Nz+(nx-1)*Nz*Ny+1+Nx*Ny*Nz,:)=0;  %Ey
            A((ny-1)*Nz+(nx-1)*Nz*Ny+1+5*Nx*Ny*Nz,:)=0;  %Hz
        end
    end
end


%{
figure; 
imagesc(abs(A)); 

figure;
imagesc(abs(Ap)); 
%}
 
end