close all; 
clear all; 
clc; 

%reads in series of transmissions and determines bandgap center as function
%of magnetic field 

p1='2';
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



shift=max(trans1);

figure;
plot(freqs1,trans1-shift,'b'); 