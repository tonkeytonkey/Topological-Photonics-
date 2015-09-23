

%close all; 
%clear all; 
%clc; 

%code takes guess point and uses it to figure out Hwant
fwant=5.57*10^9; 

fguess=5.78*10^9;  
Hs=1912;
Hguess=2160; 
gamma=2.8*10^6;

iter=20; 

maxstep=300; 
curstep=maxstep;
fcur=fguess; 
Hcur=Hguess; 

omega0=2*pi*gamma*Hguess;
omegam=2*pi*gamma*Hs;
                       
a1=1+omega0*omegam/(omega0^2-(fguess*2*pi)^2);
khr=(fguess*2*pi)*omegam/(omega0^2-(fguess*2*pi)^2);

det0=a1-khr; 

%disp(fcur);
%disp(Hcur);

for n=1:1:iter
    
    if(fwant>fcur)
        Hcur=Hcur+curstep;
    else
        Hcur=Hcur-curstep; 
    end
    curstep=curstep*3/4; 
    
    omega0=2*pi*gamma*Hcur;
    omegam=2*pi*gamma*Hs;
                       
    a1=1+omega0*omegam/(omega0^2-(fcur*2*pi)^2);
    khr=(fcur*2*pi)*omegam/(omega0^2-(fcur*2*pi)^2);

    det=a1-khr; 
    
    fcur=fguess*(1-(det-det0)/det0/2);  %this is from perturbation theory  
    
    %disp(fcur);
    %disp(round(Hcur));
    
end

disp(fcur);
disp(round(Hcur));