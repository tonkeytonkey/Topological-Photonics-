%function bandsTETMTriangle()

%close all; 
%clear all; 
%clc; 


interpolate=10; 
bands=10;  
resultnum=1; 

results='results';

nv=sprintf('%d',resultnum); 
dirsub=[results,nv];
     
if(exist(dirsub,'dir')==1)  %remove directory if already present 
    rmdir(dirsub);
end
mkdir(dirsub); 


Gamma=[0,0];
K=[0.5,0.5];
M=[0,0.5]; %we need to cover the surface BZ which is just along y axis
    
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
    
sd=size(k); 
sd=sd(1,2);

omgs=zeros(bands,sd);

rad=0.11;
a2=15;
%a1=0.84;
%khr=0.41;
a1=1;
khr=0; 

data = [];
count = 1;
for n=1:1:sd

%18,22,35,40

%n=18;  

    [omega,fem] = SquareRodTM(bands,k(1,n)*pi*2,k(2,n)*pi*2,rad,a2,a1,khr);
    %[omega,fem] = TriTM2(k(1,n)*pi*2,k(2,n)*pi*2); %need to scale k differently than we do in our original problem
    omgs(:,n)=omega;
    disp(n/sd*100); 
    data = [data;omega];
end


figure;
plot(0:1:(sd-1),data,'b-','linewidth',2);
xlim([0,(sd-1)]);
%ylim([0,1.4]);

for nk=1:1:sd

    a=sprintf('%d',nk); 
    fname=[dirsub,'/','output',a,'.txt']; 
    fp=fopen(fname,'w');
    for n=1:1:bands
        fprintf(fp,'%d \n',omgs(n,nk));
    end
    fclose(fp);
    
end



band=10;
omegamin=100;
omegamax=0; 
for n=1:1:sd
    if(omgs(band,n)>omegamax)
        omegamax=omgs(band,n);
    end
    if(omgs(band+1,n)<omegamin)
        omegamin=omgs(band+1,n);
    end
end

disp(omegamax);
disp(omegamin);
disp('gapsize');
disp(200*(omegamin-omegamax)/(omegamin+omegamax)); 


%{

kai=0.5;
%
knMG = 8*5;
knGK = 10*5;
knKM = 5*5;

MG = [];
for ii = 1:knMG
    kx = 0;
    ky = 2*pi/3^0.5-(ii-1)*2*pi/3^0.5/knMG;
    omega = TriTM(kx,ky);
    MG = [MG;omega];
end

GK = [];
for jj = 1:knGK
    kx = (jj-1)*2*pi/3/knGK;
    ky = (jj-1)*2*pi/3^0.5/knGK;
    omega = TriTM(kx,ky);
    GK = [GK;omega];
end

KM = [];
for kk = 1:knKM
    kx = 2*pi/3-(kk-1)*2*pi/3/knKM;
    ky = 2*pi/3^0.5;
    omega = TriTm(kx,ky);
    KM = [KM;omega];
end

data = [MG;GK;KM];
data = [data;data(1,:)];
figure;
%hold on
plot(data,'-ob','linewidth',2);
kntotal = knMG+knGK+knKM+1; %+1 is the last M point redundent
kTick=[1 knMG+1 knMG+knGK+1 kntotal];
set(gca,'xTick',kTick,'XTickLabel',{'M';'Gamma';'K';'M'},'FontName','Gill Sans MT','FontSize',18);
xlim([1,kntotal]);
ylim([0,1.0]);
ylabel('a/\lambda');
%tit = sprintf('',eps);
%title('Bandstructure');

%}