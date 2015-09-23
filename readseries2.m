close all; 
clear all; 
clc; 


p1={'3','4','5','6','7','9','10','11','12','13','14'};

range1=[6.17,6.09,5.99,5.99,5.9,5.55,5.55,5.55,5.55,5.55,5.53]*10^9;
range2=[6.59,6.48,6.43,6.38,6.33,6.17,6.17,6.11,5.95,5.99,5.79]*10^9;

B=[      0.305,0.295,0.290,0.281,0.275,0.260,0.259,0.256,0.248,0.250,0.239];

%these 
Btheory=[3680 ,3280 ,3050 ,2940 ,2750 ,2270 ,2270 ,2240 ,2120 ,2150 ,2020 ]/10^4;

sd=size(range1);
sd=sd(1,2); 

cent=zeros(1,sd);

day='sep5/';

for nd=1:1:sd

    tx='.txt';

    fp=fopen([day,char(p1(1,nd)),tx],'r');

    for n=1:1:2
        fgets(fp);
    end

    disp(fgets(fp)); %display comments

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
    
    
end


figure; 
plot(cent,B,'b*'); 
xlabel('frequency');
ylabel('magnetic field'); 

hold on;
plot(cent,Btheory,'r*'); 

figure;
plot(cent,Btheory-B,'g*'); 


%we want to print out B and Btheory for processing with other data sets

fp=fopen('b2.txt','w');
for nd=1:1:sd
    fprintf(fp,'%f %f \n',B(1,nd),Btheory(1,nd)); 
end
fclose(fp);






