close all; 
clear all; 
clc; 

dir='attenresults';

resultname='bob'; 
resultnum='4'; 

dir=[dir,resultnum]; 


    
    fp=fopen([dir,'/',resultname,'x','.txt'],'r');
    datax=textscan(fp,'%f %f %f ');
    fclose(fp);

    fp=fopen([dir,'/',resultname,'y','.txt'],'r');
    datay=textscan(fp,'%f %f %f ');
    fclose(fp);

    datax=cell2mat(datax); 
    datay=cell2mat(datay); 

    sd1=size(datax);
    sd1=sd1(1,1);
    posx1=zeros(1,sd1);
    posy1=zeros(1,sd1);
    valx=zeros(1,sd1);

    posx1(1,:)=datax(:,1);
    posy1(1,:)=datax(:,2);
    valx(1,:)=datax(:,3);

    sd2=size(datay);
    sd2=sd2(1,1);
    posx2=zeros(1,sd2);
    posy2=zeros(1,sd2);
    valy=zeros(1,sd2);

    posx2(1,:)=datay(:,1);
    posy2(1,:)=datay(:,2);
    valy(1,:)=datay(:,3);

    minx=min(posx1);
    maxx=max(posx1);

    miny=min(posy1);
    maxy=max(posy1);

    sx=1000;
    sy=round(sx*(maxy-miny)/(maxx-minx));

    x=linspace(minx,maxx,sx);
    y=linspace(miny,maxy,sy);

    px=griddata(posx1,posy1,valx,x,y');
    py=griddata(posx2,posy2,valy,x,y');  %py and px original have different data set sizes, but we interpolate them to be the same size 

    
    
    figure; 
    imagesc(x,-y,px); 
    axis equal; 
    ylim([-maxy,-miny]);
    xlim([minx,maxx]);  
    colormap(bluewhitered);
    title('px'); 

    figure; 
    imagesc(x,-y,py); 
    axis equal; 
    ylim([-maxy,-miny]);
    xlim([minx,maxx]);  
    colormap(bluewhitered);
    title('py'); 
    hold on;
    
    
%{
    
%down boundary, atten2
xmin1=3;
xmax1=9;
ymin1=4;
ymax1=27;

yrange=ymin1:1:ymax1; 
    
s=size(yrange);
sc1=s(1,2)-1; 
    
pc1=zeros(1,sc1);  %we want to calculate power as a function of distance 
    
for n=1:1:sc1
    
    xcount1=round((xmin1-minx)/(maxx-minx)*sx);
    xcount2=round((xmax1-minx)/(maxx-minx)*sx);
    ycount1=round((yrange(1,n)-miny)/(maxy-miny)*sy);
    ycount2=round((yrange(1,n+1)-miny)/(maxy-miny)*sy);

    int2=0;
    for nx=xcount1:1:xcount2
        for ny=ycount1:1:ycount2
            int2=int2+py(ny,nx);   %boundary in py direction
        end
    end
    int2=-int2/(x(1,xcount1)-x(1,xcount2))/(y(1,ycount1)-y(1,ycount2)); %going down so get minus 
    plot([x(1,xcount1),x(1,xcount1),x(1,xcount2),x(1,xcount2),x(1,xcount1)],-[y(1,ycount1),y(1,ycount2),y(1,ycount2),y(1,ycount1),y(1,ycount1)],'g');
    pc1(1,n)=int2; 
    
end
   
c=polyfit(yrange(1,1:sc1),log(pc1),1); 

yrange2=ymin1:(0.1):ymax1;
s=size(yrange2);
sfit=s(1,2); 
pfit=yrange2*c(1,1)+c(1,2); 
pfit=exp(pfit); 

figure;
plot(yrange(1,1:sc1),pc1,'b*'); 
hold on; 
plot(yrange2,pfit,'g'); 
    
disp('decay parameter in a');
disp(1/c(1,1)); 

%110 periods 
   
 %}   
    
    %{
    
%left boundary, atten3
xmin1=5;
xmax1=27;
ymin1=4;
ymax1=8;

xrange=xmin1:1:xmax1; 
    
s=size(xrange);
sc1=s(1,2)-1; 
    
pc1=zeros(1,sc1);  %we want to calculate power as a function of distance 
    
for n=1:1:sc1
    
    xcount1=round((xrange(1,n)-minx)/(maxx-minx)*sx);
    xcount2=round((xrange(1,n+1)-minx)/(maxx-minx)*sx);
    ycount1=round((ymin1-miny)/(maxy-miny)*sy);
    ycount2=round((ymax1-miny)/(maxy-miny)*sy);

    int2=0;
    for nx=xcount1:1:xcount2
        for ny=ycount1:1:ycount2
            int2=int2-px(ny,nx);   %boundary in px direction
        end
    end
    int2=-int2/(x(1,xcount1)-x(1,xcount2))/(y(1,ycount1)-y(1,ycount2)); %going down so get minus 
    plot([x(1,xcount1),x(1,xcount1),x(1,xcount2),x(1,xcount2),x(1,xcount1)],-[y(1,ycount1),y(1,ycount2),y(1,ycount2),y(1,ycount1),y(1,ycount1)],'g');
    pc1(1,n)=int2; 
    
end
   
c=polyfit(xrange(1,1:sc1),log(pc1),1); 

xrange2=xmin1:(0.1):xmax1;
s=size(xrange2);
sfit=s(1,2); 
pfit=xrange2*c(1,1)+c(1,2); 
pfit=exp(pfit); 

figure;
plot(xrange(1,1:sc1),pc1,'b*'); 
hold on; 
plot(xrange2,pfit,'g'); 
    
disp('decay parameter per a');
disp(1/c(1,1)); 

%190

%}

    
    
%left boundary, atten4
xmin1=14;
xmax1=38;
ymin1=4;
ymax1=8;

xrange=xmin1:1:xmax1; 
    
s=size(xrange);
sc1=s(1,2)-1; 
    
pc1=zeros(1,sc1);  %we want to calculate power as a function of distance 
    
for n=1:1:sc1
    
    xcount1=round((xrange(1,n)-minx)/(maxx-minx)*sx);
    xcount2=round((xrange(1,n+1)-minx)/(maxx-minx)*sx);
    ycount1=round((ymin1-miny)/(maxy-miny)*sy);
    ycount2=round((ymax1-miny)/(maxy-miny)*sy);

    int2=0;
    for nx=xcount1:1:xcount2
        for ny=ycount1:1:ycount2
            int2=int2-px(ny,nx);   %boundary in px direction
        end
    end
    int2=-int2/(x(1,xcount1)-x(1,xcount2))/(y(1,ycount1)-y(1,ycount2)); %going down so get minus 
    plot([x(1,xcount1),x(1,xcount1),x(1,xcount2),x(1,xcount2),x(1,xcount1)],-[y(1,ycount1),y(1,ycount2),y(1,ycount2),y(1,ycount1),y(1,ycount1)],'g');
    pc1(1,n)=int2; 
    
end
   
c=polyfit(xrange(1,1:sc1),log(pc1),1); 

xrange2=xmin1:(0.1):xmax1;
s=size(xrange2);
sfit=s(1,2); 
pfit=xrange2*c(1,1)+c(1,2); 
pfit=exp(pfit); 

figure;
plot(xrange(1,1:sc1),pc1,'b*'); 
hold on; 
plot(xrange2,pfit,'g'); 
    
disp('decay parameter per a');
disp(1/c(1,1)); 

%360


