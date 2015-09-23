


%we feed in parameters from the computational side of things which just
%include the number of processors and the directory to read the data into

close all; 
clear all; 
clc; 

    dirname='results';
    dirnum=3;
    interpolate=20;
    
    
    dirsub=[dirname,sprintf('%d',dirnum),'/'];
    if(exist(dirsub,'dir')==0)
        mkdir(dirsub)
    end
    
    
    %specify the k-points we want to investigate in the system

    M=[1/2,0];
    Gamma=[0,0];    
    K=[1/2,1/2]; 

    %interpolate=8; 
    kpoints=4; 
    list=zeros(2,kpoints); 
    list(:,1)=Gamma;
    list(:,2)=M;
    list(:,3)=K;
    list(:,4)=Gamma; 

    k=zeros(2,kpoints+(kpoints-1)*interpolate); 

    for n=1:1:(kpoints-1)
        for n2=1:1:2
            k(n2,((n-1)*(interpolate+2)+1-1*(n-1)):(n*(interpolate+2)-1*(n-1)))=linspace(list(n2,n),list(n2,n+1),interpolate+2); 
        end
    end
    
    sk=size(k);
    sk=sk(1,2); 
    
    
    %{
    dispbands=[4,5];
    freqref=(0.8421+0.7964)/2; 
    rad=0.15; 
    %}
    
    %{
    dispbands=[5,6];
    freqref=(1.109+1.071)/2; 
    rad=0.105; 
    %}
    
    dispbands=[6,7];
    freqref=(1.1664+1.1361)/2; 
    rad=0.1025; 
    
    freqscale=14.5*10^9;  %we place feature where we want inbetween the two bands 
    maxint=20; 
   
    sb=size(dispbands);
    sb=sb(1,2); 
    
    bands=dispbands(1,sb); %we only need to get bands up to this one
            
    for nband=1:1:sb
   
        fprintf('band %d \n',nband); 
        bandnum=sprintf('%d',nband);
        
        for nk=1:1:sk 
        
            fprintf('\n \n \n \n \n');
            
            knum=sprintf('%d',nk);
            fpout=fopen([dirsub,'/output',knum,',',bandnum,'.txt'],'w');
            
            fprintf('k %d \n',nk);
            freqprev=freqref; %we approximate the frequency as the frequency between the two bands we're interested in to start, then we plug in successively better guesses back in
             
            skip=-1; 
            
            for nint=1:1:maxint
        
                
                %we need to create a new settings file every iteration, we
                %don't need any text outputs from the solver 
                
                if(skip==-1)
                
                    freqstore=freqprev;
                    
                    fprintf('iteration %d \n',nint); 

                    %these correspond to experimentally realizable parameters 
                    a2=15;
            
                    f=freqprev/freqref*freqscale;  %results output from code are already frequencies 
                    gamma=2.8*10^6;
                    Hs=1780;
                    H0=2020;
            
                    omega0=2*pi*gamma*H0;
                    omegam=2*pi*gamma*Hs;
            
            
                    a1=1+omega0*omegam/(omega0^2-(f*2*pi)^2);
                    khr=(f*2*pi)*omegam/(omega0^2-(f*2*pi)^2); 
            
                    fprintf('mu11 and mu12*I %f %f ',a1,khr);
        
                    [omegas,fem]=SquareRodTM(bands,k(1,nk)*pi*2,k(2,nk)*pi*2,rad,a2,a1,khr);
        
                    omegas=sort(real(omegas));
                
                    freqprev=omegas(dispbands(1,nband)); 
                
                else
                    
                    freqprev=freqstore; 
                    fprintf('skip'); 
                    
                end
                
                fprintf('frequency %f \n',freqprev);
                fprintf(fpout,'%f \n',freqprev); %we record the frequences as we do subsequent iterations
                
                %if results are converged to one part in 10^6, we stop
                %iterating, this will reduce run time for many cases 
                if(abs((freqprev-freqstore)/freqstore)<10^-7)
                    skip=1; 
                end
            
            end
            
            fclose(fpout);
            
        end
    end
    
    

    