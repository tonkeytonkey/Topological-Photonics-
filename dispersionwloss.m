close all; 
clear all; 
clc; 


%we want to plot the loss uxx and uxy to make sure it makes physical sense 



factor=0.5:0.001:1.5; 

s=size(factor);
s=s(1,2); 

a1real=zeros(1,s);
a1imag=zeros(1,s); 

for n=1:1:s

    gamma=2.8*10^6;
    Hs=1780;
    H0=2020;
    dH=40;  %this is from Pozar  
   
    w0=2*pi*gamma*H0;
    wm=2*pi*gamma*Hs;
    w=w0*factor(1,n); 
    alpha=dH*gamma/2/w; 

    w0=w0+i*alpha*w; 
    a1=1+w0*wm/(w0^2-w^2);
    khr=w*wm/(w0^2-w^2);

    a1real(1,n)=real(a1); 
    a1imag(1,n)=imag(a1); 
    
end

figure; 
plot(factor,a1real);

figure; 
plot(factor,a1imag); 