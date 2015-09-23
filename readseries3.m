close all; 
clear all; 
clc; 


p1={'4e','5e','6e','7e','8e','9e','10e','11e','12e','13e','14e','15e','16e','17e','18e','20e'};

range1=[4.28,4.28,4.41,4.41,4.49,4.49,4.49,4.49,4.67,4.67,4.54,4.54,4.69,4.69,4.74,4.74]*10^9;
range2=[4.78,4.89,4.92,4.92,4.95,4.97,4.99,4.99,5.14,5.02,5.02,4.99,5.02,5.02,5.02,5.02]*10^9;

B=[     0.204,0.216,0.227,0.231,0.238,0.243,0.247,0.250,0.260,0.260,0.250,0.252,0.269,0.269,0.283,0.287];

%these 
Btheory=[1690 ,1820 ,2000 ,2010 ,2180 ,2270,2320 ,2330, 3400, 2800, 2400, 2400, 2800, 2800, 3090, 3090   ]/10^4;

sd=size(range1);
sd=sd(1,2); 

cent=zeros(1,sd);

day='sep9/';

for nd=1:1:sd

    tx='.txt';

    fp=fopen([day,char(p1(1,nd)),tx],'r');

    for n=1:1:2
        fgets(fp);
    end

    disp(fgets(fp)); %display comments
    disp(B(1,nd)); 
    
    for n=1:1:5
        fgets(fp);  
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqs1=data(:,1)';
    trans1=data(:,2)';
    sd1=size(freqs1);
    sd1=sd1(1,1); 

    fclose(fp); 
    
    [val,ind1]=min(abs(freqs1-range1(1,nd)));
    [val,ind2]=min(abs(freqs1-range2(1,nd)));
    
    wsum=0;
    sum=0;
    for n=ind1:1:ind2
        sum=sum+trans1(1,n);
        wsum=wsum+trans1(1,n)*freqs1(1,n); 
    end
    
    cent(1,nd)=wsum/sum; 
    
    shift=max(trans1);

    figure;
    plot(freqs1,trans1-shift,'b'); 
    hold on; 
    
    plot(freqs1(1,ind1:ind2),trans1(1,ind1:ind2)-shift,'g');
    plot([cent(1,nd),cent(1,nd)],[min(trans1),max(trans1)]-shift,'r'); 
    
    gamma=2.8*10^6;        
    f0=gamma*Btheory(1,nd)*10^4; %convert to gauss
    plot([f0,f0],[min(trans1),max(trans1)]-shift,'k'); 
   
    xlim([min(freqs1),max(freqs1)]); 
    
end



figure; 
plot(cent,B,'b*'); 
xlabel('frequency');
ylabel('magnetic field'); 



hold on;
plot(cent,Btheory,'r*'); 

freqs=min(cent):(0.01*max(cent)):max(cent); 

plot(freqs,freqs/gamma/10^4,'k*'); 

xlim([min(cent),max(cent)]); 


fp=fopen('b3.txt','w');
for nd=1:1:sd
    fprintf(fp,'%f %f \n',B(1,nd),Btheory(1,nd)); 
end
fclose(fp);

