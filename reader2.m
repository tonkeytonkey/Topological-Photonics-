close all; 
clear all; 
clc; 

%read in processed data so can run faster 

dir='results'; 
resultnum='2';
dir=[dir,resultnum]; 

fp=fopen(['output',resultnum,'.txt'],'r');
data=textscan(fp,'%f %f %f %f \n');
fclose(fp);
data=cell2mat(data);
data=data';
s=size(data);
len=s(1,2);

hm=zeros(1,len);
t=zeros(1,len);
t1=zeros(1,len);
t2=zeros(1,len);

hm(1,:)=data(1,:);
t(1,:)=data(2,:);
t1(1,:)=data(3,:);
t2(1,:)=data(4,:);

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