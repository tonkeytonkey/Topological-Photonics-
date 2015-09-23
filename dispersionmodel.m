close all;
clear all;
clc; 


freqspace=1000; 

freqs=linspace(0.5,13.5,freqspace)*10^9; 

k=zeros(1,freqspace); 
mu=zeros(1,freqspace); 

for n=1:1:freqspace

    freq=freqs(1,n); 
    
    dH=17; %resonance width
    losstan=0.0002; 
    gamma=2.8*10^6;
        
    Hs=1275;
    H0=0; 
           
    omega0=2*pi*gamma*H0;
    omegam=2*pi*gamma*Hs;
    omega0=omega0+i*(2*pi*gamma*dH)/2; 
        
    mu(1,n)=1+omega0*omegam/(omega0^2-(freq*2*pi)^2);
    k(1,n)=(freq*2*pi)*omegam/(omega0^2-(freq*2*pi)^2); 
    
end

figure;
plot(freqs,real(k)); 


figure;
plot(freqs,imag(k)); 