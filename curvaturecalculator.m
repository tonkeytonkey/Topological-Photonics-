

%we want to read in the epsilon files and critical tensor files and make
%sure everything is working properly and reassemble the bandstructure
%output from the split up k points

close all; 
clear all; 
clc; 

 
%{
runnumber=1;
name='epssweep';
sweeprange=9; 
subdir='sweep';
%}

%{
runnumber=1;
name='rodsweep';
sweeprange=5; 
subdir='rodsweep';
%}


runnumber=3;
name='rodsweep';
sweeprange=1; 
subdir='sweep';

sort=1; 

%{
runnumber=3;
name='insuldesign';
sweeprange=1; 
subdir='sweep';
%}

%{
runnumber=1;
name='design2/rodsweep';
sweeprange=1; 
subdir='sweep';
%}

%we need to isolate a chunk of the BZ containing a single dirac cone 

gridsize=10; 
    
kx=linspace(0,1,gridsize); 
ky=linspace(0,1,gridsize); 
    
k=zeros(2,gridsize*gridsize); 

for n1=1:1:(gridsize)
    for n2=1:1:(gridsize)
        k(1,(n1-1)*gridsize+n2)=kx(1,n1); 
        k(2,(n1-1)*gridsize+n2)=ky(1,n2); 
    end
end

s=size(k);
sk=s(1,2);

%determine the sizes of all the buffers

a=sprintf('%d',runnumber); 

%we need to determine the size of the system for readin

fp=fopen([name,a,'/',subdir,'1/pos.txt'],'r');
sx=fscanf(fp,'%d',1);
sy=fscanf(fp,'%d',1); 
disp('size of system is sx by sy ');
disp(sx);
disp(sy);
x=zeros(1,sx); 
for n=1:1:sx
    x(1,n)=fscanf(fp,'%f',1); 
end
y=zeros(1,sy); 
for n=1:1:sy
    y(1,n)=fscanf(fp,'%f',1); 
end
fclose(fp); 
    
dk=kx(1,2)-kx(1,1);

bands=1; 

offset=0; 

phase1=zeros(1,gridsize-1,gridsize-1);
phase2=zeros(1,gridsize-1,gridsize-1);
eigenvectors=zeros(2,2,3*sx*sy); 

%we need to read in the B matrix which gives us our inner product

B=speye(3*sx*sy,3*sx*sy)*0; 
fp=fopen([name,a,'/',subdir,'1/Breal.txt'],'r'); 
si=fscanf(fp,'%d ',1);  
for n1=1:1:si
    out=fscanf(fp,'%f %f %f ',3); 
    id=out(1,1); 
    jd=out(2,1); 
    val=out(3,1); 
    B(id,jd)=val+B(id,jd); 
end
fclose(fp); 

fp=fopen([name,a,'/',subdir,'1/Bimag.txt'],'r'); 
si=fscanf(fp,'%d ',1); 
for n1=1:1:si
    out=fscanf(fp,'%f %f %f ',3); 
    %disp(out); 
    id=out(1,1); 
    jd=out(2,1); 
    val=out(3,1); 
    B(id,jd)=1i*val+B(id,jd); 
end
fclose(fp); 


%we assume the inner product is working, in other versions we validated
%this, we assume the eigenvectors are already sorted so we can just read
%the one we want from the output file directly 

sumtot=0; 

datareal=zeros(1,3*sx*sy*bands);
dataimag=zeros(1,3*sx*sy*bands); 

for nd=4:1:4
    
    phasetot=0; 
    
    %add up little phase loops at each point in brillouin zone and sum to
    %find total phase loops
    for nkx=1:1:(gridsize-1)
        for nky=1:1:(gridsize-1)
            
            for na=1:1:2
                for nb=1:1:2
    
                    %pick correct kvector file
                    
                    kind=(nky-2+nb)*gridsize+(nkx+na-1);
                    b=sprintf('%d',kind); 
                    
                    %disp(kind); 
                    
                    %next pick correct level
                    if(mod(k(1,kind),1)==0 && mod(k(2,kind),1)==0)
                        if(nd==1)
                            disp('hi'); 
                            eigenvectors(na,nb,:)=ones(1,3*sx*sy); %we need to stick in a dummy eigenvector 
                        else
                            c=sprintf('%d',nd-1); 
                                %read in eigenvector file
                            fp=fopen([name,a,'/',subdir,'1/eigenvectorreal',b,'vec',c,'.txt'],'r');
                            datareal=fscanf(fp,'%f ',3*sx*sy); 
                            fclose(fp); 
                    
                            fp=fopen([name,a,'/',subdir,'1/eigenvectorimag',b,'vec',c,'.txt'],'r');
                            dataimag=fscanf(fp,'%f ',3*sx*sy); 
                            fclose(fp); 
                            
                            eigenvectors(na,nb,:)=datareal(:)+1i*dataimag(:);
                        end
                        
                    else
                        
                        c=sprintf('%d',nd); 
                                %read in eigenvector file
                        fp=fopen([name,a,'/',subdir,'1/eigenvectorreal',b,'vec',c,'.txt'],'r');
                        datareal=fscanf(fp,'%f ',3*sx*sy); 
                        fclose(fp); 
                    
                        fp=fopen([name,a,'/',subdir,'1/eigenvectorimag',b,'vec',c,'.txt'],'r');
                        dataimag=fscanf(fp,'%f ',3*sx*sy); 
                        fclose(fp); 
                    
                        eigenvectors(na,nb,:)=datareal(:)+1i*dataimag(:);
                    end
                    
                    
                end
            end
            
            %return; 
            
            p1=reshape(eigenvectors(1,1,:),3*sx*sy,1)'*B*reshape(eigenvectors(2,1,:),3*sx*sy,1);
            p1=p1/abs(p1); 
            p2=reshape(eigenvectors(2,1,:),3*sx*sy,1)'*B*reshape(eigenvectors(2,2,:),3*sx*sy,1);
            p2=p2/abs(p2); 
            p3=reshape(eigenvectors(2,2,:),3*sx*sy,1)'*B*reshape(eigenvectors(1,2,:),3*sx*sy,1);
            p3=p3/abs(p3); 
            p4=reshape(eigenvectors(1,2,:),3*sx*sy,1)'*B*reshape(eigenvectors(1,1,:),3*sx*sy,1);
            p4=p4/abs(p4); 
            
            %loops that include dirac cones should yield plus or minus pi
            %phase 
            
            prod=imag(log(p2))-imag(log(p4'))-imag(log(p3'))+imag(log(p1));
            
            phase2(1,nkx,nky)=prod; 
            
            prod=p1*p2*p3*p4;
            prod=imag(log(prod)); 

            phasetot=phasetot+prod; 
            phase1(1,nkx,nky)=prod;
            
            
           
        end
    end
    
    sumtot=phasetot/(2*pi)+sumtot; 
    
    disp('band number');
    disp(nd);
    disp(phasetot/(2*pi));
    disp(sumtot); 

end
    


for nd=1:1:1
    figure;
    imagesc(reshape(phase1(nd,:,:),gridsize-1,(gridsize)-1)); 
    axis equal; 
    %surf(kx,ky,reshape(omegas(nd,:,:),gridsize,gridsize));
    a=sprintf('%d',nd);
    title(a); 
end
for nd=1:1:1
    figure;
    imagesc(reshape(phase2(nd,:,:),gridsize-1,(gridsize)-1)); 
    axis equal; 
    %surf(kx,ky,reshape(omegas(nd,:,:),gridsize,gridsize));
    a=sprintf('%d',nd);
    title(a); 
end

sum=0;
for nx=1:1:3
    for ny=3:1:6
        sum=phase2(1,nx,ny)+sum;
    end
end

disp('circulation around one dirac cone'); 
disp(sum);

sum=0;
for nx=3:1:6
    for ny=1:1:3
        sum=phase2(1,nx,ny)+sum;
    end
end

disp('circulation around one dirac cone'); 
disp(sum);

sum=0;
for nx=3:1:6
    for ny=6:1:9
        sum=phase2(1,nx,ny)+sum;
    end
end

disp('circulation around one dirac cone'); 
disp(sum);

sum=0;
for nx=6:1:9
    for ny=3:1:6
        sum=phase2(1,nx,ny)+sum;
    end
end

disp('circulation around one dirac cone'); 
disp(sum);

%we want to write little script to check if bands touching

