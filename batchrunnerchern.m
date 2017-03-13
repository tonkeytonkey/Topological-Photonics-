


%we feed in parameters from the computational side of things which just
%include the number of processors and the directory to read the data into

function [] = batchrunnerchern(dirmain,split,splitnum)

%{
    close all; 
    clear all; 
    clc; 
    dirmain='testresults'; 
    split=1;
    splitnum=1; 
%}
    
    I=i; 

    sx=64; 
    sy=64; 
    sz=1; 
    omegamin=0.01; 
    omegamax=1.6;   %we need to keep all the bands we're interested in 
    center=1.1*(2*pi); 
    
    bands=40; 
    
    %subname='rodsweep'; 
    subname='sweep';
    
    
    %specify the k-points we want to investigate in the system

    gridsize=10; 
    
    kx=linspace(0,1,gridsize); 
    kx=kx(1,1:gridsize); %we don't want to include boundary cells twice
    ky=linspace(0,1,gridsize); 
    ky=ky(1,1:gridsize);
    
    k=zeros(2,gridsize*gridsize); 

    for n1=1:1:gridsize
        for n2=1:1:gridsize
            k(1,(n1-1)*gridsize+n2)=kx(1,n1); 
            k(2,(n1-1)*gridsize+n2)=ky(1,n2); 
        end
    end

    kind=1:1:(gridsize*gridsize);  %we need to save eigenvectors at all points
    
    svary=1; 
    
    for nvary=1:1:svary
        
        disp('param step'); 
        disp(nvary);
        
        nv=sprintf('%d',nvary); 
        dirsub=[dirmain,'/',subname,nv];
      
        %we first check if the directory exists, and if it doesn't we make
        %it for this particular parameter run
        if(exist(dirsub,'dir')==0)
            mkdir(dirsub)
        end
      
    %now specify the system geometry
    
        R1=[1,0,0];
        R2=[0,1,0]; 
        R3=[0,0,1]; 

        Repeatx=1; 
        Repeaty=1; 
        Repeatz=1; 

        Mx=0; 
        My=0; 
        Mz=0; 
        
        symblen=2; 

        symbolcoord=zeros(symblen,3); 
        symbolvec=zeros(symblen,3); 
        symbolrep=zeros(symblen,1);
        symboltype=zeros(symblen,1); 
        symbolparam=zeros(symblen,2); 
        symbolten=zeros(symblen,4,3,3); 
    
        %we want to create a plain square lattice of rods

         
        
        symbolcount=1; 

        symbolcoord(symbolcount,:)=[0,0,0];     %gave in cartesian coordinates, not lattice coordinates, if repeats more than once, the repeat vector connects to the neighbors
        symbolvec(symbolcount,:)=0; %give the repeat vector
        symbolrep(symbolcount,:)=1; %gives the number of times we repeat the symbol 
        symboltype(symbolcount,1)=2; %two means block 
        symbolparam(symbolcount,1)=norm(R1);
        symbolparam(symbolcount,2)=norm(R2);
        symbolten(symbolcount,1,:,:)=[1,0,0; 0,1,0; 0,0,1]; %delta,epsilon
        symbolten(symbolcount,2,:,:)=[1,0,0; 0,1,0; 0,0,1];    %delta,m
        symbolten(symbolcount,3,:,:)=0;    %other guys
        symbolten(symbolcount,4,:,:)=0; 
        symbolcount=symbolcount+1; 

         %a1=eps(1,nvary); 
        a2=14.28; 

        f=6.6*10^9; 
        
        gamma=2.8*10^6;
        
        Hs=1817;
        H0=2400; 
           
        omega0=2*pi*gamma*H0;
        omegam=2*pi*gamma*Hs;
        
        a1=1+omega0*omegam/(omega0^2-(f*2*pi)^2);
        khr=(f*2*pi)*omegam/(omega0^2-(f*2*pi)^2); 
        
        r=0.324/3.00; 
        
        
        symbolcoord(symbolcount,:)=R1/2+R2/2;     %gave in cartesian coordinates, not lattice coordinates, if repeats more than once, the repeat vector connects to the neighbors
        symbolvec(symbolcount,:)=0; %give the repeat vector
        symbolrep(symbolcount,:)=1; %gives the number of times we repeat the symbol 
        symboltype(symbolcount,1)=1; %one means circle
        symbolparam(symbolcount,1)=r; %set both sidelengths of the square
        symbolten(symbolcount,1,:,:)=[a2,0,0; 0,a2,0; 0,0,a2]; %delta,epsilon
        symbolten(symbolcount,2,:,:)=[a1,khr*I,0; -khr*I,a1,0; 0,0,1];    %delta,mu
        symbolten(symbolcount,3,:,:)=0;    %other guys
        symbolten(symbolcount,4,:,:)=0; 
        symbolcount=symbolcount+1; 

     
        
        s=size(symbolparam); 
        sparam=s(1,2); 

        
        
        %we create one geometry file for each individual thread, this
        %prevents weird read conflicts
        
        %if(exist([dirsub,'/settings.txt'],'file')==0)
        
         d=sprintf('%d',splitnum); 
        
         fp=fopen([dirsub,'/settings',d,'.txt'],'w');  %print out position information
         %first need to write out position vectors and repeat distances
         fprintf(fp,'lattice vectors \n'); 
         fprintf(fp, '   %d %d %d   \n',R1(1,1),R1(1,2),R1(1,3)); 
         fprintf(fp, '   %d %d %d   \n',R2(1,1),R2(1,2),R2(1,3));
         fprintf(fp, '   %d %d %d   \n',R3(1,1),R3(1,2),R3(1,3));
         fprintf(fp,'repeat numbers \n'); 
         fprintf(fp,'%d \n',Repeatx); 
         fprintf(fp,'%d \n',Repeaty); 
         fprintf(fp,'%d \n',Repeatz); 
         fprintf(fp,'metal boundary \n'); 
         fprintf(fp,'%d \n',Mx); 
         fprintf(fp,'%d \n',My); 
         fprintf(fp,'%d \n',Mz); 
         fprintf(fp,'symbollen \n %d \n',symbolcount-1); 
         fprintf(fp,'symbfields \n %d \n',sparam); 
    
         fprintf(fp,'type, params,   coord   ,copy vec, rep  ,  tensors \n');
            
         for n=1:1:(symbolcount-1)
            fprintf(fp,' %d    ',symboltype(n,1)); 
            for n2=1:1:sparam
                fprintf(fp, '%d ',symbolparam(n,n2)); %we may have multiple symbol parameters
            end
            fprintf(fp, '   %d %d %d    ',symbolcoord(n,1),symbolcoord(n,2),symbolcoord(n,3)); 
            fprintf(fp, '   %d %d %d    ',symbolvec(n,1),symbolvec(n,2),symbolvec(n,3));
            fprintf(fp, '%d   ',symbolrep(n,1)); 
            for n2=1:1:4
                for na=1:1:3
                    for nb=1:1:3
                        fprintf(fp,'%d ',real(symbolten(n,n2,na,nb))); %output all tensor information for the symbol
                    end
                 end
                 fprintf('   '); 
             end
             fprintf(fp,'\n                                 '); 
             for n2=1:1:4
                for na=1:1:3
                    for nb=1:1:3
                        fprintf(fp,'%d ',imag(symbolten(n,n2,na,nb))); %output all tensor information for the symbol
                    end
                end
                fprintf('   '); 
             end
                
             fprintf(fp,' \n'); %go to the next line
         end
        fclose(fp); 
      

        %now that we've created the file we can actually run the job
        
        %bianistropicsolver1(dirsub,sx,sy,sz,k,bands,center,omegamin,omegamax,kind,split,splitnum); 
        bianistropicsolver1TE(dirsub,sx,sy,sz,k,bands,center,omegamin,omegamax,kind,split,splitnum); 


        if(splitnum~=1)
            delete([dirsub,'/settings',d,'.txt']); %we want to clean up the directory after we're done
        end
        
        
        
    end
    
    
return; 
    