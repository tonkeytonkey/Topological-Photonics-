close all; 
clear all; 
clc; 

%purpose of this code is to read in the mode profiles 

dir='results3';

resultname='bob'; 
len=2; 

t=zeros(1,len);
t1=zeros(1,len); 
t2=zeros(1,len); 
hm=zeros(1,len); 

for n=1:1:len

    b=sprintf('%d',n); 
    
    fp=fopen([dir,'/',resultname,b,'E','.txt'],'r');
    name=fgets(fp); %we can get rid of annoying line at the top, this contains hm
    hm(1,n)=cell2mat(textscan(name,'%f'));
    dataE=textscan(fp,'%f %f %f ');
    fclose(fp);

    dataE=cell2mat(dataE); 

    sd1=size(dataE);
    sd1=sd1(1,1);
    posx1=zeros(1,sd1);
    posy1=zeros(1,sd1);
    valE=zeros(1,sd1);

    posx1(1,:)=dataE(:,1);
    posy1(1,:)=dataE(:,2);
    valE(1,:)=dataE(:,3);

    minx=min(posx1);
    maxx=max(posx1);

    miny=min(posy1);
    maxy=max(posy1);

    sx=1000;
    sy=round(sx*(maxy-miny)/(maxx-minx));

    x=linspace(minx,maxx,sx);
    y=linspace(miny,maxy,sy);

    Ez=griddata(posx1,posy1,valE,x,y');

    %we want to plot the circles 
    points=100;
    H=15;
    W=28; 
    offset=0.20; 
    offset2=0.1; 
    shift=0; 
    a=0.805; 
    wm=0.5; 
    Wpml=2; 
    Hpml=1.8;
    
    %{
    m1=rect2(wm,hm,'base','corner','pos',{'5.75',H-hm},'rot','0');


%creates left lattice
g4=circ2('r1','base','center','pos',{'0.5',0.5+offset2},'rot','0','const',fem.const);
gt1=geomarrayr(g4,1,1,W/2+shift,H);

%creates right lattice 
g5=circ2('r2','base','center','pos',{0.5+W/2+shift,0.5+offset},'rot','0','const',fem.const);
gt2=geomarrayr(g5,a,a,round((W/2-shift)/a),round(H/a)-1);
%}
    figure; 
    hold on; 
    
    imagesc(x,y,Ez); 
    
    colormap(bluewhitered);
    
    xpoints=zeros(1,points);
    ypoints=zeros(1,points); 
    
    r=0.15; 
    xcen=0.5;
    ycen=0.5+offset2; 
    for nx=1:1:(W/2+shift)
        for ny=1:1:H
            for np=1:1:points
                xpoints(1,np)=r*cos((np-1)/points*2*pi)+xcen;
                ypoints(1,np)=r*sin((np-1)/points*2*pi)+ycen;
            end
            xpoints(1,points+1)=xpoints(1,1);
            ypoints(1,points+1)=ypoints(1,1); 
            plot(xpoints,ypoints,'k','LineWidth',1);
        ycen=ycen+1; 
        end
        xcen=xcen+1; 
        ycen=0.5+offset2; 
    end
    
    xcen=0.5+W/2+shift;
    ycen=0.5+offset; 
    for nx=1:1:round((W/2-shift)/a)
        for ny=1:1:(round(H/a)-1)
            for np=1:1:points
                xpoints(1,np)=r*cos((np-1)/points*2*pi)+xcen;
                ypoints(1,np)=r*sin((np-1)/points*2*pi)+ycen;
            end
            xpoints(1,points+1)=xpoints(1,1);
            ypoints(1,points+1)=ypoints(1,1); 
            plot(xpoints,ypoints,'k','LineWidth',1);
        ycen=ycen+a; 
        end
        xcen=xcen+a; 
        ycen=0.5+offset; 
    end
    
    cx=5.75;
    cy=H-hm(1,n); 
    hn=0.1; 
    
    fill([0,cx,cx,cx+wm,cx+wm,W-Wpml-4,W-Wpml-4,0,0],[H,cy+hm(1,n),cy,cy,cy+hm(1,n),H,H+hn,H+hn,H],[0.7,0.7,0.7],'LineWidth',1); 
    
    axis equal; 
    ylim([Hpml+5.6,H+hn]); 
    xlim([0,W-Wpml-4.0]); 
    set(gca,'xTick',[]); 
    set(gca,'yTick',[]); 

    annotation('textbox',...
    [0.3 0.43 0.075 0.075],...
    'String',{'C=1'},...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','-',...
    'EdgeColor',[1 1 1],...
    'LineWidth',2,...
    'BackgroundColor',[1  1 1],...
    'Color',[0 0 0]);
    
    annotation('textbox',...
    [0.75 0.43 0.075 0.075],...
    'String',{'C=2'},...
    'FontSize',14,...
    'FontName','Arial',...
    'LineStyle','-',...
    'EdgeColor',[1 1 1],...
    'LineWidth',2,...
    'BackgroundColor',[1  1 1],...
    'Color',[0 0 0]);
    
    %ylim([miny,maxy]);
    %xlim([minx,maxx]);  
    
end




%{
figure;
hold on;
plot(hm,t,'g*-');
plot(hm,t1,'b*-');
plot(hm,t2,'r*-'); 
ylim([0,1]); 
legend('Total','Channel 1','Channel 2'); 

fz=20;
lw=0.5; 

set(get(gca,'Ylabel'),'String','Power (a.u.)','FontSize',fz,'FontName','Arial','LineWidth',lw); 
set(get(gca,'Xlabel'),'String','Metal Height (a.u.)','FontSize',fz,'FontName','Arial','LineWidth',lw);
set(get(gca,'Xlabel'),'FontSize',fz,'FontName','Arial','LineWidth',lw); 
set(get(gca,'Ylabel'),'FontSize',fz,'FontName','Arial','LineWidth',lw); 
set(gca,'FontSize',fz,'FontName','Arial','LineWidth',lw);
set(gca,'xTick',[1, 1.25, 1.5,1.75,2.0]); 
set(gca,'yTick',[0, 0.5, 1.0]); 
%}





