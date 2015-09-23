close all; 
clear all; 
clc; 

%we want to read in data sets from Cgap=1 crystal and Cgap=2 crystal to
%show that they are well aligned 

smoothwin=3; 


%feb13 for edge transmission data Cgap=2
%dec29, Cgap=1

%feb6 for Cgap=2
%dec30 for Cgap=1

ncomment=ones(1,8)*1; 

%best data sets


%read={'4S24.txt','4S42.txt'};
%names={'S24','S42'};



%read={'4S14.txt','4S41.txt'};
%names={'S14','S41'};



%0.39 T
read={'feb13\73S12.txt','dec29\73S12.txt'};
names={'C_{gap}=1','C_{gap}=2'};
xticks=[7.0,7.5,8.0,8.5];
xlim1=7.2;
xlim2=8.9; 


%{
%0.35 T
read={'feb13\65S12.txt','dec29\65S12.txt'};
names={'C_{gap}=1','C_{gap}=2'};
xticks=[7.0,7.5,8.0,8.5];
xlim1=6.9;
xlim2=8.7; 
%}

%{
%0.315 T
read={'feb13\59S12.txt','dec29\58S12.txt'};
names={'C_{gap}=1','C_{gap}=2'};
xticks=[7.0,7.5,8.0,8.5];
xlim1=6.7;
xlim2=8.5; 
%}

color={'r','b'};



%{
ncomment=1;
read={'1S14.txt'}; 
%}

%{
names={'a','b','c','d','e','f','g'}; 
color={'r','g','b','c','k','y','k'}; 
%}

s=size(read);
s=s(1,2); 

figure;
hold on; 

for n=1:1:s


    fp=fopen(read{1,n},'r');

    for n2=1:1:2
        fgets(fp);
    end

    for n2=1:1:ncomment(1,n)
        str=fgets(fp);
        disp(str); 
    end
    
    for n2=1:1:5
        fgets(fp);  
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqs1=data(:,1)';
    trans1=data(:,2)';
    
    fclose(fp); 
    
    disp(color{1,n}); 
    plot(freqs1/10^9,smooth(trans1,smoothwin)',color{1,n});
    ylabel('Transmission (dB)');
    xlabel('Frequency (GHz)'); 
    
end
    
fz=10;
lw=0.5; 

set(gca,'FontSize',fz,'FontName','Arial','LineWidth',lw);
set(gca,'xtick',xticks);
set(gca,'ytick',[-100,-70,-40]);


legend(names);
ylim([-105,-35]); 
xlim([xlim1,xlim2]); 

    
%daspect([1,120,1]); 
    
box on; 
 
    
 




