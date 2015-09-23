close all; 
clear all; 
clc; 

smoothwin=3; 

%this data is for small rods in small lattice using antenna #2 

%{
readone='jan21/large1S12.txt';
readtwo='jan21/smallbulkS12.txt'; 
readthree='jan21/large1S12.txt'; 
%}

%{

%large crystal, isolated S12

readone='jan21/large1S12.txt';
one='0.33 T large crystal'; 
readtwo='jan21/large2S12.txt'; 
two='0.35 T large crystal'; 
readthree='jan21/large3S12.txt'; 
three='0.40 T large crystal'; 
%}

%{

%small crystal isolated S12

readone='jan21/small3S12.txt';
one='0.33 T small crystal'; 
readtwo='jan21/smallbulkS12.txt'; 
two='0.35 T small crystal'; 
readthree='jan21/smallbulk2S12.txt'; 
three='0.40 T small crystal'; 
%}

%{
%comparision of 0.35 T measurement on Jan21

readone='jan21/large2S12.txt'; 
one='0.35 T large crystal'; 
readtwo='jan21/smallbulkS12.txt'; 
two='0.35 T small crystal'; 
readthree='jan21/smallbulk2S12.txt'; 
three='0.40 T small crystal'; 
%}

%{
%comparison of old and new for small crystal 
readone='jan16/5.txt';
one='0.35 T small crystal old'; 
readtwo='jan21/smallbulkS12.txt'; 
two='0.35 T small crystal new '; 
readthree='jan21/smallbulk2S12.txt'; 
three='0.40 T small crystal'; 
%}

readone='jan21/2S13.txt';
one='S13'; 
readtwo='jan21/2S31.txt'; 
two='S31 '; 
readthree='jan21/smallbulk2S12.txt'; 
three='0.40 T small crystal'; 



%{
readone='jan16/9.txt';
readtwo='dec22/76S12.txt'; 
readthree='jan16/8.txt'; 
%}

%{
readone='jan16/8.txt';
readtwo='dec29/75S12.txt'; 
readthree='jan16/8.txt'; 
%}

%{
readone='jan16/2S12.txt';
readtwo='jan16/2S21.txt'; 
readthree='jan16/8.txt'; 
%}

    fp=fopen(readone,'r');

    for n=1:1:2
        fgets(fp);
    end

    str=fgets(fp);
    disp(str); 
    
    for n=1:1:5
        fgets(fp);  
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqs1=data(:,1)';
    trans1=data(:,2)';
    
    fclose(fp); 
    
 
    
    
    
    fp=fopen(readtwo,'r');

    for n=1:1:2
        fgets(fp);
    end

    str=fgets(fp);
    disp(str); 
    
    for n=1:1:5
        fgets(fp);  
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqs2=data(:,1)';
    trans2=data(:,2)';
    
    fclose(fp); 
    
    
    
    
    fp=fopen(readthree,'r');

    for n=1:1:2
        fgets(fp);
    end

    str=fgets(fp);
    disp(str); 
    
    for n=1:1:5
        fgets(fp);  
    end

    data=textscan(fp,'%f, %f  \n');
    data=cell2mat(data); 

    freqs3=data(:,1)';
    trans3=data(:,2)';
    
    fclose(fp); 
    
    
    figure;
    hold on;
    plot(freqs1,smooth(trans1,smoothwin)','r');
    plot(freqs2,smooth(trans2,smoothwin)','b'); 
    %plot(freqs3,smooth(trans3,smoothwin)','g'); 
    legend(one,two,three); 
    ylim([-120,0]); 
    xlim([2*10^9,10^10]); 




