close all; 
clear all; 
clc; 

readstart=1;
readfinish=22; %give file numbers we want to read in
points=201; 
day='oct1/';
letter=''; 
tx='.txt';

nd=-1;
avg=-1000; 


for nd=readstart:1:readfinish
    
    textfile=[day,sprintf('%d',nd),letter,tx];
    fp=fopen(textfile,'r');

    for n=1:1:2
        fgets(fp);
    end

    disp(nd);
    disp(fgets(fp)); 
    
    for n=1:1:5
        fgets(fp);  
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqs1=data(:,1)';
    trans1=data(:,2)';

    fclose(fp); 
    
    figure;
    hold on;
    plot(freqs1,trans1,'b'); 
    title(sprintf('%d',nd));
        
    xlim([min(freqs1),max(freqs1)]); 
    
    %we want to find the best average from 3 to 10 GHz
    
    [val,ind1]=min(abs(3*10^9-freqs1));
    [val,ind2]=min(abs(10*10^9-freqs1));
    
    tot=0;
    for n=ind1:1:ind2
        tot=tot+trans1(1,n);
    end
    
    tot=tot/(ind2-ind1+1); 
    
    plot(freqs1(ind1:1:ind2),tot,'g'); 
    
    disp(tot);     
    
end





