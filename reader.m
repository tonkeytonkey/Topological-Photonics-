close all; 
clear all; 
clc; 

dir='results';

resultname='bob'; 
resultnum='1';
len=41; 

dir=[dir,resultnum]; 

t=zeros(1,len);
t1=zeros(1,len); 
t2=zeros(1,len); 
hm=zeros(1,len); 

for n=1:1:len

    b=sprintf('%d',n); 
    
    fp=fopen([dir,'/',resultname,b,'x','.txt'],'r');
    name=fgets(fp); %we can get rid of annoying line at the top, this contains hm
    hm(1,n)=cell2mat(textscan(name,'%f'));
    datax=textscan(fp,'%f %f %f ');
    fclose(fp);

    fp=fopen([dir,'/',resultname,b,'y','.txt'],'r');
    fgets(fp); 
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

    %{
    
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
    
    %}

    %need to calculate power in different slices and normalize by size
   
    %input boundary 
    xmin1=9;
    xmax1=12;
    ymin1=10;
    ymax1=15;

    xcount1=round((xmin1-minx)/(maxx-minx)*sx);
    xcount2=round((xmax1-minx)/(maxx-minx)*sx);
    ycount1=round((ymin1-miny)/(maxy-miny)*sy);
    ycount2=round((ymax1-miny)/(maxy-miny)*sy);

    int1=0;
    for nx=xcount1:1:xcount2
        for ny=ycount1:1:ycount2
            int1=int1+px(ny,nx);  %boundary in px direction
        end
    end
    int1=int1/(x(1,xcount1)-x(1,xcount2))/(y(1,ycount1)-y(1,ycount2)); 
    plot([x(1,xcount1),x(1,xcount1),x(1,xcount2),x(1,xcount2),x(1,xcount1)],-[y(1,ycount1),y(1,ycount2),y(1,ycount2),y(1,ycount1),y(1,ycount1)],'g');

    %down boundary 
     xmin1=12;
    xmax1=17;
    ymin1=5;
    ymax1=10;

    xcount1=round((xmin1-minx)/(maxx-minx)*sx);
    xcount2=round((xmax1-minx)/(maxx-minx)*sx);
    ycount1=round((ymin1-miny)/(maxy-miny)*sy);
    ycount2=round((ymax1-miny)/(maxy-miny)*sy);

    int2=0;
    for nx=xcount1:1:xcount2
        for ny=ycount1:1:ycount2
            int2=int2+py(ny,nx);   %boundary in py direction
        end
    end
    int2=-int2/(x(1,xcount1)-x(1,xcount2))/(y(1,ycount1)-y(1,ycount2)); %going down so get minus 
    plot([x(1,xcount1),x(1,xcount1),x(1,xcount2),x(1,xcount2),x(1,xcount1)],-[y(1,ycount1),y(1,ycount2),y(1,ycount2),y(1,ycount1),y(1,ycount1)],'g');

    %right boundary 
    xmin1=17;
    xmax1=26;
    ymin1=10;
    ymax1=15;

    xcount1=round((xmin1-minx)/(maxx-minx)*sx);
    xcount2=round((xmax1-minx)/(maxx-minx)*sx);
    ycount1=round((ymin1-miny)/(maxy-miny)*sy);
    ycount2=round((ymax1-miny)/(maxy-miny)*sy);

    int3=0;
    for nx=xcount1:1:xcount2
        for ny=ycount1:1:ycount2
            int3=int3+px(ny,nx);   %boundary in px direction
        end
    end
    int3=int3/(x(1,xcount1)-x(1,xcount2))/(y(1,ycount1)-y(1,ycount2)); 
    plot([x(1,xcount1),x(1,xcount1),x(1,xcount2),x(1,xcount2),x(1,xcount1)],-[y(1,ycount1),y(1,ycount2),y(1,ycount2),y(1,ycount1),y(1,ycount1)],'g');

    disp(int1); 
    disp(int3+int2);

    t(1,n)=(int2+int3)/int1;
    t1(1,n)=int2/int1;
    t2(1,n)=int3/int1; 
    
end


fp=fopen(['output',resultnum,'.txt'],'w');
for n=1:1:len
    fprintf(fp,'%f %f %f %f \n',hm(1,n),t(1,n),t1(1,n),t2(1,n));
end
fclose(fp);

figure;
hold on;
plot(hm,t,'g-','LineWidth',2);
plot(hm,t1,'b-','LineWidth',2);
plot(hm,t2,'r-','LineWidth',2); 
ylim([0,1]); 
legend('Total','Channel 1','Channel 2'); 

fz=20;
lw=0.5; 

set(get(gca,'Ylabel'),'String','Transmission','FontSize',fz,'FontName','Arial','LineWidth',lw); 
set(get(gca,'Xlabel'),'String','Metal Height (a)','FontSize',fz,'FontName','Arial','LineWidth',lw);
set(get(gca,'Xlabel'),'FontSize',fz,'FontName','Arial','LineWidth',lw); 
set(get(gca,'Ylabel'),'FontSize',fz,'FontName','Arial','LineWidth',lw); 
set(gca,'FontSize',fz,'FontName','Arial','LineWidth',lw);
set(gca,'xTick',[1, 1.25, 1.5,1.75,2.0]); 
set(gca,'yTick',[0, 0.5, 1.0]); 

%}


%{

%now we want to calculate the Hx and Hy fields, and finally the power 
Hx=circshift(Ez,[0,1])-circshift(Ez,[0,-1]);
Hy=circshift(Ez,[1,0])-circshift(Ez,[-1,0]);

figure; 
imagesc(gridx,-gridy,Hx); 
axis equal; 
ylim([-maxy,-miny]);
xlim([minx,maxx]); 
colormap(bluewhitered);

figure; 
imagesc(gridx,-gridy,Hy); 
axis equal; 
ylim([-maxy,-miny]);
xlim([minx,maxx]); 
colormap(bluewhitered);

power=((Hx.*Ez).^2+(Hy.*Ez).^2).^(1/2); 

figure; 
imagesc(gridx,-gridy,power); 
axis equal; 
ylim([-maxy,-miny]);
xlim([minx,maxx]); 
%colormap(bluewhitered);

%}







