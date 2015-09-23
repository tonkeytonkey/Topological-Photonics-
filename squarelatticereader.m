

%we want to read in the epsilon files and critical tensor files and make
%sure everything is working properly and reassemble the bandstructure
%output from the split up k points

close all; 
clear all; 
clc; 

interpolate=40; 
 
%{
runnumber=1;
name='epssweep';
sweeprange=9; 
subdir='sweep';
%}

%{
runnumber=1;
name='rodsweep';
sweeprange=5; 
subdir='rodsweep';
%}

runnumber=1;
name='results';

%{

P1=[0,-0.5];
P2=[0,0.5]; %we need to cover the surface BZ which is just along y axis
    
%interpolate=8; 
kpoints=2; 
list=zeros(2,kpoints); 
list(:,1)=P1;
list(:,2)=P2;

%}


M=[0.5,0.5]; 
Gamma=[0,0]; 
X=[0.5,0]; 

%interpolate=8; 
kpoints=4; 
list=zeros(2,kpoints); 
list(:,1)=Gamma;
list(:,2)=X;
list(:,3)=M;
list(:,4)=Gamma; 

k=zeros(2,kpoints+(kpoints-1)*interpolate); 

for n=1:1:(kpoints-1)
    for n2=1:1:2
        k(n2,((n-1)*(interpolate+2)+1-1*(n-1)):(n*(interpolate+2)-1*(n-1)))=linspace(list(n2,n),list(n2,n+1),interpolate+2); 
    end
end

s=size(k);
sk=s(1,2);

%determine the sizes of all the buffers

a=sprintf('%d',runnumber); 

nsweep=1;

    b=sprintf('%d',nsweep); 
    
    fname=[name,a,'/output1.txt'];
    disp(fname); 
    fp=fopen(fname,'r');
    data=textscan(fp,'%d');
    data=cell2mat(data);
    s=size(data);
    bands=s(1,1);
    disp('number of bands reading in ');
    disp(bands);
    fclose(fp);


    omegas=zeros(bands,sk); 

    for n=1:1:sk
        c=sprintf('%d',n); 
        fp=fopen([name,a,'/output',c,'.txt'],'r');
        data=textscan(fp,'%f \n');
        data=cell2mat(data);
        omegas(:,n)=sort(data(:,1));
        fclose(fp); 
    end

    
    bandsk=8; 
    
    keepbands=zeros(bandsk,sk); 
    
    for n=1:1:bandsk
        for nk=1:1:sk
            keepbands(n,nk)=omegas(n,nk);
        end
    end
    
    gy=[0.5,0.5,0.5];
    
    figure; 
    hold on;
    plot(0:1:(sk-1),keepbands,'k','LineWidth',2); 
    plot(0:1:(sk-1),keepbands(6,:),'Color',[0,0,1],'LineWidth',4); 
    plot(0:1:(sk-1),keepbands(7,:),'Color',[0,0,1],'LineWidth',4); 
    xlim([0,sk-1]); 
    line([interpolate+1,interpolate+1],[0,1.2],'Color',gy,'LineWidth',1);
    line([2*interpolate+2,2*interpolate+2],[0,1.2],'Color',gy,'LineWidth',1);

    fz=14;
    lw=0.5; 
    
    hmin=0;
    hmax=1.2;
    
    %xlabel('K vector'); 
    %ylim([0.4,0.8]); 
    movegui(gcf,'northwest'); 
    %ylim([0,1.2]); 
    ylim([hmin,hmax]);
    
    
    set(gca,'xtick',[]);
    set(get(gca,'Ylabel'),'String','Frequency (2\pi{}c/a)','FontSize',fz,'FontName','Arial','LineWidth',lw); 
    set(gca,'FontSize',fz,'FontName','Arial','LineWidth',lw);
    
    ax=axis; 
    
    toff=-0.25; 
    thei=-0.025+hmin; 
    
    text(toff,thei,'\Gamma','FontSize',fz,'FontName','Arial');
    text(interpolate+1+toff*2,thei,'M','FontSize',fz,'FontName','Arial');
    text(2*interpolate+2+toff*2,thei,'K','FontSize',fz,'FontName','Arial');
    text(3*interpolate+3+toff,thei,'\Gamma','FontSize',fz,'FontName','Arial');
    
    %set(gca,'YTick',[0.7 0.8 0.9 1.0 ]);
    
    set(gca,'plotboxaspectratio',[1,2,1]);
    
    box on; 
    
    %we also want to add chern numbers directly to lines

    
    
    %set(gca,'plotboxaspectratio',[1 2 1])
    
    %now we want min and max values of 2nd and 3rd bands
    
    omegamax=-1; 
    omegamin=1000;
    
    omega1=zeros(1,sk); 
    omega2=zeros(1,sk); 
    
    band1=4;
    band2=5;
    
    for nk=1:1:sk
        if(omegamax<omegas(band1,nk))
            omegamax=omegas(band1,nk); 
        end
        if(omegamin>omegas(band2,nk))
            omegamin=omegas(band2,nk); 
        end
        omega1(1,nk)=omegas(band1,nk); 
        omega2(1,nk)=omegas(band2,nk); 
    end
    
    disp('min and max omega of 3rd and 4th bands'); 
    disp(omegamin); 
    disp(omegamax);
    
    
    
    
    rectangle('Position',[0,omegamax,sk,(omegamin-omegamax)],'FaceColor','y','EdgeColor','y');
    plot(0:1:(sk-1),omegas(6,:),'b','LineWidth',4); 
    plot(0:1:(sk-1),omegas(7,:),'b','LineWidth',4); 
    
    group=[1,1,1,1,2]; 
    chern=[0,1,-1,-1,4];

    bandstart=cumsum(group); %gives number where we start our calculations from
    bandstart=circshift(bandstart,[0,1]);
    bandstart=bandstart+1;
    bandstart(1,1)=1; %gives number we start counting at for each group
    
    s=size(group);
    sg=s(1,2);
    for n=1:1:sg
        ch=sprintf('%d',chern(1,n));
        avg=0;    
        for nb=1:1:group(1,n)
            avg=avg+keepbands(nb-1+bandstart(1,n),interpolate/2+1);
        end
        avg=avg/group(1,n); 
    
        text(interpolate/2,avg,ch,'FontSize',fz,'FontName','Arial','backgroundcolor',[1,1,1]); 
    end
    
    
    
    %set(gca,'plotboxaspectratio',[1 2 1])
    
    %printeps(gcf, 'squarelattice2'); 
    
    
    disp('gapsize'); 
    disp(2*(omegamin-omegamax)/(omegamax+omegamin)*100); 
    
    figure;
    plot(1:1:sk,omega1,'r'); 
    hold on; 
    plot(1:1:sk,omega2,'b'); 
    
    
    
   








