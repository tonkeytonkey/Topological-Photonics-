close all; 
clear all; 
clc; 


p1={'16a','17a','18a','19a','20a','21a','22a','23a','24a','25a','26a','27a','28a','29a','31a','32a','33a','34a'};

range1=[4.45,4.45,4.45,4.37,4.37,4.37,4.37,4.29,4.29,4.14,4.02,3.91,3.79,3.67,3.72,3.695,3.415,3.58]*10^9;
range2=[4.93,4.93,4.93,4.87,4.84,4.84,4.81,4.71,4.68,4.68,4.74,4.61,4.54,4.11,4.00,3.955,3.645,3.82]*10^9;

B=      [0.233,0.231,0.228,0.223,0.220,0.214,0.208,0.201,0.198,0.190,0.186,0.179,0.168,0.161,0.154,0.151,0.139,0.147];

%these 
Btheory=[2050 ,2050 ,2050 ,1900 ,1860 ,1830 ,1790 ,1630 ,1610 ,1510 ,1460 ,1350 ,1250 ,1080 ,1050 ,1030 ,850 ,950 ]/10^4;

sd=size(range1);
sd=sd(1,2); 

smoothwin=1; 
cent=zeros(1,sd);

%we want to store interpolated versions of the data sets for direct
%comparison to each other 
freqref=(3:0.05:6.0)*10^9;
Bref=min(B):0.005:max(B); 

%we store each frequency and B field associated with each point 
freqstore=zeros(1,201*sd); 
Bstore=zeros(1,201*sd); 
transtore=zeros(1,201*sd); 

day='sep5/'; 

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
    trans1=smooth(trans1,smoothwin)'; 
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
    title(sprintf('%f',B(1,nd)));
    
    gamma=2.8*10^6;        
    f0=gamma*Btheory(1,nd)*10^4; %convert to gauss
    plot([f0,f0],[min(trans1),max(trans1)]-shift,'k'); 
   
    xlim([min(freqs1),max(freqs1)]); 
    
    transtore(1,(1+201*(nd-1)):(201*nd))=-(trans1-shift)/min(trans1-shift); 
    freqstore(1,(1+201*(nd-1)):(201*nd))=freqs1; 
    Bstore(1,(1+201*(nd-1)):(201*nd))=B(1,nd); 
    
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

%now we interpolate our scattered data set so we can plot using contour
%plot or other similiar functions 
[x,y]=meshgrid(freqref,Bref);
intep=griddata(freqstore,Bstore,transtore,x,y); 

figure;
imagesc(freqref,Bref,intep); 
hold on; 
plot(cent,B,'k*');
plot(Btheory*gamma*10^4,B,'w*');

set(gca, 'CLim', [-1, 0]);


%{
figure;
%imagesc(freqref,Btheory,transref);
contour(freqref,Btheory,transref); 
%}




fp=fopen('b1.txt','w');
for nd=1:1:sd
    fprintf(fp,'%f %f \n',B(1,nd),Btheory(1,nd)); 
end
fclose(fp);





