



close all; 
clear all; 
clc; 

 


runnumber=2;
name='mpbresult';

a=sprintf('%d',runnumber); 
dir=[name,a,'/']; 

gridsize=200; 
bands=10; 
    
kx=linspace(0,1,gridsize); 
ky=linspace(0,1,gridsize); 
    
k=zeros(2,gridsize*gridsize); 

for n1=1:1:gridsize
    for n2=1:1:gridsize
        k(1,(n1-1)*gridsize+n2)=kx(1,n1); 
        k(2,(n1-1)*gridsize+n2)=ky(1,n2); 
    end
end

s=size(k);
sk=s(1,2);

dk=kx(1,2)-kx(1,1); 

%first we need to check that we can read in eigenvectors and that we can
%take inner product 

%{

num1=1;
num2=2;

if(num1/10 < 1)
    b1=sprintf('0%d',num1);
else
    b1=sprintf('%d',num1);
end
if(num2/10 < 1)
    b2=sprintf('0%d',num2);
else
    b2=sprintf('%d',num2);
end

fname=['trianglesweep-h.k',b1,'.b',b2,'.z.te.h5'];

eig1=hdf5read([dir,fname],'z.r')+i*hdf5read([dir,fname],'z.i');

val=reshape(eig1,[],1)'*reshape(eig1,[],1); 

s=size(eig1);
seig=s(1,1)^2; 

%}



%we also read in bands so we can make sure they look how we expect 

bandstore=10; 

omegas=zeros(bandstore,gridsize,gridsize); 

fname=[name,a,'/','band1.dat'];
data=dlmread(fname,',',1,1); 
vals=data(:,6:end); 
kxs=data(:,2);
kys=data(:,3); 


s=size(vals);
sw=s(1,1); 

for n=1:1:sw
   omegas(:,round(kxs(n,1)/dk+1),round(kys(n,1)/dk+1))=vals(n,1:bandstore); 
end





%we want to write little script to check if bands touching

gaps=ones(1,bands-1)*100;
for nd=1:1:(bands-1)
    for nx=1:1:gridsize
        for ny=1:1:gridsize
            if((omegas(nd+1,nx,ny)-omegas(nd,nx,ny))<gaps(1,nd))
                gaps(1,nd)=(omegas(nd+1,nx,ny)-omegas(nd,nx,ny));
            end
        end
    end
end

disp('minimum gaps between bands');
disp(gaps); 




%break into two sets

%clear vars eigenvectors; 

%{
for nd=1:1:bands
    a=sprintf('%d',nd);
    figure; 
    imagesc(reshape(phase(nd,:,:),gridsize-1,gridsize-1)); 
    title(a);
end
%}


hmin=0.68;
hmax=0.95;

figure;
hold on; 
grid on;

h1=surf(kx,ky,reshape(omegas(4,:,:),gridsize,gridsize));
h2=surf(kx,ky,reshape(omegas(5,:,:),gridsize,gridsize));
shading interp;

h3=pcolor(kx,ky,reshape(log(0.005+abs(omegas(4,:,:)-omegas(5,:,:))),gridsize,gridsize));
shading interp;
set(h3,'ZData',hmin+0*reshape(omegas(4,:,:)-omegas(5,:,:),gridsize,gridsize)); %have to put this weird statement in so the size of Zdata doesn't change

[C,h4]=contour(kx,ky,reshape(omegas(4,:,:)-omegas(5,:,:),gridsize,gridsize),0,'LineColor',[1,1,1]*0);

plot3([0,0.5,0,0],[0,0.5,0.5,0],[hmin,hmin,hmin,hmin],'--','Color',[1,1,1]*1,'LineWidth',3);

%even though caxis is the same for both objects, we can modify the cdata
%for each object so that we end up with appropriate scaling 

newmap=[linspace(1,0,64);linspace(0,0,64);linspace(0,1,64)].';

colormap([flipud(jet(256));jet(128);fliplr(jet(128))]);
%colormap([flipud(jet(256));jet(256)]);

c1=get(h1,'CData');
c2=get(h2,'CData');
c1new=(c1-0.7)/(hmax-0.7)/2+0.5;   %map from [0.5,1] hsv range of colormap
c2new=(c2-0.7)/(hmax-0.7)/2+0.5;
set(h1,'CData',c1new);
set(h2,'CData',c2new);

c3=get(h3,'CData'); 
c3min=min(min(c3));
c3max=max(max(c3));
c3new=(c3-c3min)/(c3max-c3min)/2;  %map from [0,0.5] gray range of colormap
set(h3,'CData',c3new);

caxis([0,1]);
%alpha(0.9); 

zlim([hmin,0.9]); 

%need to change fontsize, type, linewidth

fz=36;
lw=0.5; 

set(get(gca,'Xlabel'),'FontSize',fz,'FontName','Arial','LineWidth',lw); 
set(get(gca,'Ylabel'),'FontSize',fz,'FontName','Arial','LineWidth',lw); 
set(get(gca,'Zlabel'),'FontSize',fz,'FontName','Arial','LineWidth',lw); 
set(gca,'FontSize',fz,'FontName','Arial','LineWidth',lw);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

offset=0.02;
%text(0,0,hmin-offset,'\Gamma','FontSize',fz,'FontName','Arial');
%text(0,0.5,hmin-offset,'M','FontSize',fz,'FontName','Arial');
%text(0.5,0,hmin-offset,'M','FontSize',fz,'FontName','Arial');

%material dull; 

set(gca,'zTick',[0.7,0.8 0.9 1.0]); 


box on; 

view([-16,12]);
%camlight(-82,34); 
%camlight(60,40);
