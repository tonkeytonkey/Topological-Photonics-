close all; 
clear all; 
clc; 


p1={'21f','22f','23f','24f','25f','26f','27f','28f','29f','30f','31f','32f','33f','34f','35f','36f','37f','38f','39f','40f','41f','42f','43f','44f','45f','46f','47f','48f','49f','50f','51f','52f','53f'};

range1=[5.96,5.96,5.96,6.03,6.03,6.243,6.27,6.27,6.325,6.38,6.54,6.54,6.10,6.01,6.01,6.01,5.88,5.88,5.73,5.83,5.83,5.83,5.76,5.76,5.73,5.68,5.47,5.45,5.375,5.375,5.45,5.56,5.77]*10^9;
range2=[6.57,6.57,6.63,6.63,6.63,6.573,6.6,6.6,6.8,    6.9, 6.98,6.98,6.44,6.42,6.42,6.42,6.38,6.34,6.05,6.05,6.19,6.19,5.96,6.06,5.9 ,5.88,5.66,5.66,5.450,5.450,5.54,5.74,6.03]*10^9;

B=      [0.287,0.291,0.293,0.295,0.297,0.308,0.310,0.312,0.319,0.330,0.337,0.340,0.294,0.289,0.289,0.286,0.275,0.271,0.258,0.257,0.267,0.267,0.255,0.260,0.249,0.247,0.239,0.235,0.227,0.228,0.231,0.241,0.256];

%these 
Btheory=[3190 ,3250 ,3350 ,3450 ,3500 ,3760 ,3970 ,4010 ,4800 ,5500 ,5500 ,5500 ,3250 ,3070 , 3070, 3070, 2800, 2730, 2330, 2400, 2500, 2500, 2280 ,2360 ,2200 ,2160 ,1920 ,1900 ,1770 ,1770 ,1850, 2010 ,2340]/10^4;

sd=size(range1);
sd=sd(1,2); 

cent=zeros(1,sd);

smoothwin=1; 

%we want to store interpolated versions of the data sets for direct
%comparison to each other 
freqref=(5:0.05:7.0)*10^9;
Bref=min(B):0.005:max(B); 

%we store each frequency and B field associated with each point 
freqstore=zeros(1,201*sd); 
Bstore=zeros(1,201*sd); 
transtore=zeros(1,201*sd); 

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
%contourf(freqref,Bref,intep,5); 
imagesc(freqref,Bref,intep);
hold on; 
plot(cent,B,'k*');
plot(Btheory*gamma*10^4,B,'w*');

set(gca, 'CLim', [-1, 0]);

%we want to print out B and Btheory for processing with other data sets

fp=fopen('b4.txt','w');
for nd=1:1:sd
    fprintf(fp,'%f %f \n',B(1,nd),Btheory(1,nd)); 
end
fclose(fp);






