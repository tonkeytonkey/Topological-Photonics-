close all; 
clear all; 
clc; 

%we give pairs of files we want to compare 

%additional data from 35b to 38b for surface transmission

%{
%C=-1 gap, 0.197 T,crappy decreased distance of metal waveguide
p1='39b';
p2='40b'; 
%}

%{
%C=-1 gap, 0.197 T, very deep minimum (all after increased distance of metal waveguide away from crystal)
p1='41b';
p2='42b'; 
%}

%{
%C=-1 gap, 0.197 T, side minimum
%interesting 
p1='43b';
p2='44b'; 
%}

%{
%C=-1 gap, 0.197 T, best example, 40 to 50 dB different 
p1='45b';
p2='46b'; 
%}


%C=2 gap, 0.274 T, best example, over 30 dB different 
p1='47b';
p2='48b'; 


%{
%C=2 gap, 0.274 T
p1='49b';
p2='50b'; 
%}

tx='.txt';

fp=fopen([p1,tx],'r');

for n=1:1:2
    fgets(fp);
end

disp(fgets(fp)); %display comments

for n=1:1:5
   fgets(fp);  
end

data=textscan(fp,'%f, %f  \n');
data=cell2mat(data); 

freqs1=data(:,1);
trans1=data(:,2);
sd=size(freqs1);
sd1=sd(1,1); 

fclose(fp); 

fp=fopen([p2,tx],'r');

for n=1:1:2
    fgets(fp);
end

disp(fgets(fp)); %display comments

for n=1:1:5
   fgets(fp);  
end

data=textscan(fp,'%f, %f  \n');
data=cell2mat(data); 

freqs2=data(:,1);
trans2=data(:,2);
sd=size(freqs2);
sd2=sd(1,1); 

fclose(fp);

shift=max([trans1;trans2]);


figure;
plot(freqs1,trans1-shift,'b'); 

hold on; 

plot(freqs2,trans2-shift,'r'); 



