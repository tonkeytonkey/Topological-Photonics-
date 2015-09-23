close all; 
clear all; 
clc; 

fp=fopen('boundary.txt','r');

data=textscan(fp,'%f %f %f');
data=cell2mat(data);

sd=size(data);
sd=sd(1,1);

x=zeros(1,sd);
y=zeros(1,sd);
val=zeros(1,sd);

x=data(:,1);
y=data(:,2);
val=data(:,3); 

minx=min(x);
maxx=max(x);

miny=min(y);
maxy=max(y);

sx=300;
sy=100;

gridx=linspace(minx,maxx,sx);
gridy=linspace(miny,maxy,sy);

gridval=griddata(x,-y,val,gridx,gridy');

h1=figure; 
imagesc(gridx,gridy,gridval); 
hold on; 

xcen=0; 
ycen=-6.8; 
r=0.07;
points=50;
xpoints=zeros(1,points+1); %need +1 so we can connect it back to itself
ypoints=zeros(1,points+1);

%want to plot array of circles 
for n2=1:1:15
    
    for n=1:1:points
        xpoints(1,n)=r*cos((n-1)/points*2*pi)+xcen;
        ypoints(1,n)=r*sin((n-1)/points*2*pi)+ycen;
    end
    xpoints(1,points+1)=xpoints(1,1);
    ypoints(1,points+1)=ypoints(1,1); 
    plot(xpoints,ypoints,'k');

    ycen=ycen+1; 
    
end

%need to also draw circles 

colormap(bluewhitered);
vmax=max(gridval(:));
vmin=min(gridval(:)); 

set(gca,'XTickLabel','');
set(gca,'YTickLabel','');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

axis equal;
%xlim([-0.5,0.5]);
%ylim([miny,0]); 

%{
ti = get(gca,'TightInset'); 
set(gca,'Position',[ti(1) ti(2) 1-ti(3)-ti(1) 1-ti(4)-ti(2)]);

set(gca,'units','centimeters')
pos = get(gca,'Position');
ti = get(gca,'TightInset');

set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize', [pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition',[0 0 pos(3)+ti(1)+ti(3) pos(4)+ti(2)+ti(4)]);
%}

%axis tight; 
colorbar; 

cmax=max(-vmin,vmax);
caxis([vmin,vmax]);  %make absolute maximum value the same on both sides of the scale 

fclose(fp);

printeps(h1,'one'); 

%print(h1,'-djpeg','-r600','one.png'); 