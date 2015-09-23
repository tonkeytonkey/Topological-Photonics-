close all; 
clear all; 
clc; 


interpolate=20; 
bandnum=350;  
resultnum=10; 

results='results';

nv=sprintf('%d',resultnum); 
dirsub=[results,nv];
     
if(exist(dirsub,'dir')==1)  %remove directory if already present 
    rmdir(dirsub);
end
mkdir(dirsub); 


P1=[0,-0.5];
P2=[0,0];
P3=[0,0.5]; %we need to cover the surface BZ which is just along y axis
    
kpoints=3; 
list=zeros(2,kpoints); 
list(:,1)=P1;
list(:,2)=P2;
list(:,3)=P3;
    
k=zeros(2,kpoints+(kpoints-1)*interpolate); 

for n=1:1:(kpoints-1)
    for n2=1:1:2
        k(n2,((n-1)*(interpolate+2)+1-1*(n-1)):(n*(interpolate+2)-1*(n-1)))=linspace(list(n2,n),list(n2,n+1),interpolate+2); 
    end
end
    
sd=size(k); 
sd=sd(1,2);

omgs=zeros(bandnum,sd);


data = [];
count = 1;
for n=1:1:sd

%18,22,35,40

%n=18;  

    [omega,fem] = SquareRodTMWGpec(k(2,n)*pi*2);
    %[omega,fem] = SquareRodTMWG(k(2,n)*pi*2);
    omgs(:,n)=omega;
    disp(n/sd*100); 
    data = [data;omega];
end



figure;
plot(data,'b-','linewidth',2);


for nk=1:1:sd

    a=sprintf('%d',nk); 
    fname=[dirsub,'/','output',a,'.txt']; 
    fp=fopen(fname,'w');
    for n=1:1:bandnum
        fprintf(fp,'%d \n',omgs(n,nk));
    end
    fclose(fp);
    
end