

%purpose of this code is to check for bandgaps and make sure they are
%complete between neighboring bands

close all; 
clear all; 
clc; 


runnumber=3;
name='rodsweep';
sweeprange=1; 
subdir='sweep';




gridsize=40; 
    
kx=linspace(0,1,gridsize); 
ky=linspace(0,1,gridsize); 
    
k=zeros(2,gridsize*gridsize); 

for n1=1:1:gridsize
    for n2=1:1:gridsize
        k(1,(n1-1)*gridsize+n2)=kx(1,n1); 
        k(2,(n1-1)*gridsize+n2)=ky(1,n2); 
    end
end

s=size(k);
sk=s(1,2);

%determine the sizes of all the buffers

a=sprintf('%d',runnumber); 
    
dk=kx(1,2)-kx(1,1);

bands=7; 

offset=0; 

omegas=zeros(bands-offset,gridsize,gridsize); 

for n=1:1:sk
    
        b=sprintf('%d',n);
    
        %read in eigenvalue file at particular k-point 
        fp=fopen([name,a,'/',subdir,'1/eigenvalue',b,'.txt'],'r');
        
        datavalue=textscan(fp,'%f'); 
        datavalue=cell2mat(datavalue);
        omegasel=datavalue; 
        fclose(fp);
  
        
        for nd=(offset+1):1:(bands)  %we only keep first band eigenvalues
            omegas(nd-offset,round(k(1,n)/dk+1),round(k(2,n)/dk+1))=omegasel(nd,1);
        end
        
end


%{

for nd=1:1:(bands-offset)
    figure;
    imagesc(reshape(omegas(nd,:,:),gridsize,gridsize)); 
    %surf(kx,ky,reshape(omegas(nd,:,:),gridsize,gridsize));
    a=sprintf('%d',nd);
    title(a); 
end

%}

%we want to write little script to check if bands touching

gaps=ones(1,bands-1)*100;
for nd=1:1:(bands-1)
    for nx=1:1:gridsize
        for ny=1:1:gridsize
            if((omegas(nd+1,nx,ny)-omegas(nd,nx,ny))<gaps(1,nd))
                gaps(1,nd)=(omegas(nd+1,nx,ny)-omegas(nd,nx,ny));
            end
        end
    end
end

disp('minimum gaps between bands');
disp(gaps); 




%break into two sets

%clear vars eigenvectors; 



figure;
hold on; 
for nd=1:1:4
    a=sprintf('%d',nd);
    surfc(kx,ky,reshape(omegas(nd,:,:),gridsize,gridsize));
end

figure;
hold on; 
for nd=4:1:5
    a=sprintf('%d',nd);
    surfc(kx,ky,reshape(omegas(nd,:,:),gridsize,gridsize));
end
%zlim([0.6,0.9]); 

figure;
hold on; 
for nd=5:1:7
    a=sprintf('%d',nd);
    surfc(kx,ky,reshape(omegas(nd,:,:),gridsize,gridsize));
end
