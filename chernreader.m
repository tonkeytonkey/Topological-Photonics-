

%we want to read in the epsilon files and critical tensor files and make
%sure everything is working properly and reassemble the bandstructure
%output from the split up k points

close all; 
clear all; 
clc; 

 


runnumber=2;
name='rodsweep';
sweeprange=1; 
subdir='sweep';

sort=1; 



gridsize=10; 
    
kx=linspace(0,1,gridsize); 
ky=linspace(0,1,gridsize); 
    
k=zeros(2,gridsize*gridsize); 

for n1=1:1:gridsize
    for n2=1:1:gridsize
        k(1,(n2-1)*gridsize+n1)=kx(1,n1); 
        k(2,(n2-1)*gridsize+n1)=ky(1,n2); 
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



datareal=zeros(1,3*sx*sy);
dataimag=zeros(1,3*sx*sy); 

%we group together the bands that we see are touching
%groups=[1,1,1,1,1,1]; 
groups=[1,2]; 
%groups=[1,1,1,1,2]; 
s=size(groups);
sg=s(1,2); 


phase=zeros(sg,gridsize-1,gridsize-1); 
kind=zeros(4,1);
sum=0; 

bandstart=cumsum(groups); %gives number where we start our calculations from
bandstart=circshift(bandstart,[0,1]);
bandstart=bandstart+1;
bandstart(1,1)=1; %gives number we start counting at for each group



for ng=1:1:sg
    
    phasetot=0; 
    
    %add up little phase loops at each point in brillouin zone and sum to
    %find total phase loops
    for nkx=1:1:(gridsize-1)
        for nky=1:1:(gridsize-1)
            
            %we need to read in eigenvectors we need to calculating for
            %this placate 
            
            kind(1,1)=(nky-1)*gridsize+(nkx);
            kind(2,1)=(nky-1)*gridsize+(nkx+1);
            kind(3,1)=(nky)*gridsize+(nkx);
            kind(4,1)=(nky)*gridsize+(nkx+1);
            
            %kind=(nky-2+nb)*gridsize+(nkx+na-1);
            
            p1=zeros(groups(1,ng),groups(1,ng));
            p2=zeros(groups(1,ng),groups(1,ng));
            p3=zeros(groups(1,ng),groups(1,ng));
            p4=zeros(groups(1,ng),groups(1,ng));
            
            eigens=zeros(3*sx*sy,4,groups(1,ng)); 
            
            for nb=1:1:groups(1,ng)
            
                for n=1:1:4
                
                    b=sprintf('%d',kind(n,1)); 
     
                    if(mod(k(1,kind(n,1)),1)==0 && mod(k(2,kind(n,1)),1)==0)
                        if(ng==1 && nb==1) %first band of 1st group
                            %disp('hi'); 
                            eigens(:,n,nb)=ones(1,3*sx*sy); %we need to stick in a dummy eigenvector 
                            c=sprintf('%d',0);
                        else
                            c=sprintf('%d',nb+bandstart(1,ng)-2); 
                                %read in eigenvector file
                            fp=fopen([name,a,'/',subdir,'1/eigenvectorreal',b,'vec',c,'.txt'],'r');
                            datareal=fscanf(fp,'%f ',3*sx*sy); 
                            fclose(fp); 
                    
                            fp=fopen([name,a,'/',subdir,'1/eigenvectorimag',b,'vec',c,'.txt'],'r');
                            dataimag=fscanf(fp,'%f ',3*sx*sy); 
                            fclose(fp); 
                            
                            eigens(:,n,nb)=datareal(:)+1i*dataimag(:);
                        end
                        
                    else
                        
                        c=sprintf('%d',nb+bandstart(1,ng)-1); 
                                %read in eigenvector file
                        fp=fopen([name,a,'/',subdir,'1/eigenvectorreal',b,'vec',c,'.txt'],'r');
                        datareal=fscanf(fp,'%f ',3*sx*sy); 
                        fclose(fp); 
                    
                        fp=fopen([name,a,'/',subdir,'1/eigenvectorimag',b,'vec',c,'.txt'],'r');
                        dataimag=fscanf(fp,'%f ',3*sx*sy); 
                        fclose(fp); 
                    
                        eigens(:,n,nb)=datareal(:)+1i*dataimag(:);
                    end
                    
                end
                
            end
             
            %calculate link variables for groups of bands 
            for nb1=1:1:groups(1,ng)
                for nb2=1:1:groups(1,ng)
                    p1(nb1,nb2)=eigens(:,1,nb1)'*B*eigens(:,2,nb2);
                    p2(nb1,nb2)=eigens(:,2,nb1)'*B*eigens(:,4,nb2);
                    p3(nb1,nb2)=eigens(:,4,nb1)'*B*eigens(:,3,nb2);
                    p4(nb1,nb2)=eigens(:,3,nb1)'*B*eigens(:,1,nb2); 
                end
            end
                      
            p1d=det(p1)/abs(det(p1)); 
            p2d=det(p2)/abs(det(p2)); 
            p3d=det(p3)/abs(det(p3)); 
            p4d=det(p4)/abs(det(p4)); 
            
            prod=p1d*p2d*p3d*p4d;
            prod=imag(log(prod))/(2*pi); 
            %disp(prod);
            
            phasetot=phasetot+prod; 
            phase(ng,nkx,nky)=prod;
           
        end
    end
    
    sum=phasetot+sum; 
    
    disp('group number');
    disp(ng);
    disp(phasetot);
    disp(sum); 

end





    


for nd=1:1:(3)
    figure;
    imagesc(reshape(phase(nd,:,:),gridsize-1,gridsize-1)); 
    %caxis([0,1]); 
    axis equal; 
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    %caxis([-1,1]); 
    colorbar; 
    
    fz=12;
    lw=0.5; 
    
    offset=0.5;
    
    text(offset,offset,'\Gamma','FontSize',fz,'FontName','Arial');
    text((gridsize-1)/2+offset,offset,'X','FontSize',fz,'FontName','Arial');
    text(offset,(gridsize-1)/2+offset,'X','FontSize',fz,'FontName','Arial');
    text((gridsize-1)/2+offset,(gridsize-1)/2+offset,'M','FontSize',fz,'FontName','Arial');
    
    
end

%we want to write little script to check if bands touching


