close all; 
clear all; 
clc; 


%this data is for large lattice

readstart=1;
readfinish=106; %give file numbers we want to read in
points1=401; 
day='dec29/';
letter='S12';
smoothwin=5; 
calibration='none';
%calibration='../oct1/1.txt';  %this is response function for old antennas from earlier file 
freqref=(0.5:0.005:13.0)*10^9;

tx='.txt';

sd=(readfinish-readstart+1); 

B1=zeros(1,sd);
freqs1=zeros(1,points1); 
trans1=zeros(sd,points1); 

%we store each frequency and B field associated with each point 
freqstore1=zeros(1,points1*sd); 
Bstore1=zeros(1,points1*sd); 
transtore1=zeros(1,points1*sd); 

%read in reference antenna function and use to properly calibrate data 

cal=1;
if(strcmp(calibration,'none')==1)
    cal=-1; %we turn off calibration function
end

if(cal==1)

    fp=fopen(calibration,'r');

    for n=1:1:8
        fgets(fp);
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqsref=data(:,1)';
    transref=data(:,2)';
    
    fclose(fp); 

end

for nd=readstart:1:readfinish
    
    
    textfile=[day,sprintf('%d',nd),letter,tx];
    fp=fopen(textfile,'r');

    for n=1:1:2
        fgets(fp);
    end

    str=fgets(fp);
    etr=sscanf(str,'  %s %f ');
    B1(1,nd)=etr((end-1),1); 
    
    for n=1:1:5
        fgets(fp);  
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqs1=data(:,1)';
    trans1(nd,:)=data(:,2)';
    trans1(nd,:)=smooth(trans1(nd,:),smoothwin)'; 
    sd1=size(freqs1);
    sd1=sd1(1,1); 

    fclose(fp); 
    
    %find calibration function based on antenna 
    if(cal==1)
        cal1=interp1(freqsref,transref,freqs1);
    end
    
    if(cal==1)
        transtore1(1,(1+points1*(nd-1)):(points1*nd))=trans1(nd,:)-cal1;
    else
        transtore1(1,(1+points1*(nd-1)):(points1*nd))=trans1(nd,:);
    end
    freqstore1(1,(1+points1*(nd-1)):(points1*nd))=freqs1; 
    Bstore1(1,(1+points1*(nd-1)):(points1*nd))=B1(1,nd); 
    
end




%this data is for small lattice 

readstart=1;
readfinish=102; %give file numbers we want to read in
points2=401; 
day='dec22/';
letter='S12';
smoothwin=5; 
calibration='none';
%calibration='../oct1/1.txt';  %this is response function for old antennas from earlier file 
freqref=(0.5:0.005:13.0)*10^9;

tx='.txt';

sd=(readfinish-readstart+1); 

B2=zeros(1,sd);
freqs2=zeros(1,points2); 
trans2=zeros(sd,points2); 

%we store each frequency and B field associated with each point 
freqstore2=zeros(1,points2*sd); 
Bstore2=zeros(1,points2*sd); 
transtore2=zeros(1,points2*sd); 

%read in reference antenna function and use to properly calibrate data 

cal=1;
if(strcmp(calibration,'none')==1)
    cal=-1; %we turn off calibration function
end

if(cal==1)

    fp=fopen(calibration,'r');

    for n=1:1:8
        fgets(fp);
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqsref=data(:,1)';
    transref=data(:,2)';
    
    fclose(fp); 

end

for nd=readstart:1:readfinish
    
    
    textfile=[day,sprintf('%d',nd),letter,tx];
    fp=fopen(textfile,'r');

    for n=1:1:2
        fgets(fp);
    end

    str=fgets(fp);
    etr=sscanf(str,'  %s %f ');
    B2(1,nd)=etr((end-1),1); 
    
    for n=1:1:5
        fgets(fp);  
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqs2=data(:,1)';
    trans2(nd,:)=data(:,2)';
    trans2(nd,:)=smooth(trans2(nd,:),smoothwin)'; 
    sd1=size(freqs2);
    sd1=sd1(1,1); 

    fclose(fp); 
    
    %find calibration function based on antenna 
    if(cal==1)
        cal1=interp1(freqsref,transref,freqs2);
    end
    
    if(cal==1)
        transtore2(1,(1+points2*(nd-1)):(points2*nd))=trans2(nd,:)-cal1;
    else
        transtore2(1,(1+points2*(nd-1)):(points2*nd))=trans2(nd,:);
    end
    freqstore2(1,(1+points2*(nd-1)):(points2*nd))=freqs2; 
    Bstore2(1,(1+points2*(nd-1)):(points2*nd))=B2(1,nd); 
    
end

Bmin=0.30; 
Bmax=0.50; 

[val,ind1start]=min(abs(B1-Bmin));
[val,ind2start]=min(abs(B2-Bmin));

[val,ind1end]=min(abs(B1-Bmax));

%we want to plot the spectra on top of each other so we can see how well
%the bandgaps are lined up, and over what range
ind2=ind2start;
for ind1=ind1start:1:ind1end
    
    figure;
    plot(freqs1,reshape(trans1(ind1,:),[1,points1]),'b'); 
    hold on;
    plot(freqs2,reshape(trans2(ind2,:),[1,points2]),'r'); 
    title(sprintf('%f %f',B1(1,ind1),B2(1,ind2)));
    ind2=ind2+1; 
    xlim([3.0,9.0]*10^9); 
end











