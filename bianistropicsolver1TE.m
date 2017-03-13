
%computes bands for square lattice with identical set of parameters as mpb
%simulation, we control the number of interpolations, the number of bands
%and the size

%omegamin and omegamax determine the ranges of eigenstates we will save

%we want one copy of this code to ensure the integrity of all of the
%results, this will be updated at one location

function [data] = bianistropicsolver1TE(dirt,sx,sy,sz,k,bands,center,omegamin,omegamax,kind,split,splitnum)


subpixel=1; 


%{
close all; 
clear all; 
clc; 

dirt='results1';
sx=20;
sy=20; 
interpolate=5; 
bands=20; 
split=2;
splitnum=1; 
omegamin=0.15;
omegamax=0.2; 
%}

disp(dirt); 

I=i; 

d=sprintf('%d',splitnum); 

fp=fopen([dirt,'/settings',d,'.txt'],'r');  %print out position information

disp('fp handle');
disp(fp);

%need to scan in R vectors with repeatX and repeatY
R1=zeros(1,3); 
R2=zeros(1,3); 
R3=zeros(1,3); 
fgets(fp); %need to remove textline
R1(1,:)=fscanf(fp,'%f ',3); 
R2(1,:)=fscanf(fp,'%f ',3); 
R3(1,:)=fscanf(fp,'%f ',3);
fgets(fp); 
Repeatx=fscanf(fp,'%f ',1); 
Repeaty=fscanf(fp,'%f ',1); 
Repeatz=fscanf(fp,'%f ',1); 
fgets(fp); 
Mx=fscanf(fp,'%f ',1); 
My=fscanf(fp,'%f ',1); 
Mz=fscanf(fp,'%f ',1); 
fgets(fp); 
symblen=fscanf(fp,'%d');
fgets(fp); 
sparam=fscanf(fp,'%d'); 
%create matrices which store symbol information
symbolcoord=zeros(symblen,3); 
symbolvec=zeros(symblen,3); 
symbolrep=zeros(symblen,1);
symboltype=zeros(symblen,1); 
symbolparam=zeros(symblen,2); 
symbolten=zeros(symblen,4,3,3); 
fgets(fp); %read in one line which has all of the textual information
for n=1:1:(symblen)
    symboltype(n,1)=fscanf(fp,' %f ',1); 
    for n2=1:1:sparam
        symbolparam(n,n2)=fscanf(fp, '%f ',1); %we may have multiple symbol parameters
    end
    symbolcoord(n,:)=fscanf(fp, '%f %f %f ',3); 
    symbolvec(n,:)=fscanf(fp, '%f %f %f ',3);
    symbolrep(n,1)=fscanf(fp, '%d ',1); 
    for n2=1:1:4
        for na=1:1:3
            for nb=1:1:3
                symbolten(n,n2,na,nb)=fscanf(fp,'%f ',1); %output all tensor information for the symbol
            end
        end
    end
    %now we need to read in the imaginary part of the tensor
    for n2=1:1:4
        for na=1:1:3
            for nb=1:1:3
                symbolten(n,n2,na,nb)=I*fscanf(fp,'%f ',1)+symbolten(n,n2,na,nb); %output all tensor information for the symbol
            end
        end
    end
end
fclose(fp);

N1=sx*Repeatx*subpixel;
N2=sy*Repeaty*subpixel; 
N3=sz*Repeatz;


Nxyz=N1*N2*N3; 

%we can calculate the reciprocal vectors
V=abs(dot(R1,cross(R2,R3))); 

G1=2*pi*-cross(R2,R3)/V; 
G2=2*pi*cross(R3,R1)/V; 
G3=2*pi*cross(R1,R2)/V;

%now we want to create a simple geometry, create vectors of positions,
%arrays of material

r1=zeros(N1,3); 
r2=zeros(N2,3); 
a1=zeros(1,N1+1); 
a2=zeros(1,N2+1); 
%r3=zeros(N3,3); 

for n=1:1:3
    a1(1,:)=linspace(0,Repeatx*R1(1,n),N1+1); 
    a2(1,:)=linspace(0,Repeaty*R2(1,n),N2+1); 
    r1(:,n)=a1(1:N1); 
    r2(:,n)=a2(1:N2); 
    
    %r3(:,n)=linspace(0,R3(1,n),N3); 
end

dx=r1(2,1)-r1(1,1);
dy=r2(2,2)-r2(1,2);


epsilonexp=zeros(N1,N2,3,3); 
muexp=zeros(N1,N2,3,3);
xiexp=zeros(N1,N2,3,3); 
zetaexp=zeros(N1,N2,3,3); 

add=zeros(N1,N2,3,3,4); 
geom=zeros(N1,N2); 
gradx=zeros(N1,N2,2,2); 
grady=zeros(N1,N2,2,2); 

%this step should establish the geometry of the system

%we need to fill in shapes being careful to fill in shape by going around
%to other side of boundary

%we want to rethink the way we are rendering object

tic; 

 for nsym=1:1:symblen
  
       if(symboltype(nsym,1)==2)  %this makes blocks 
           
            w=symbolparam(nsym,1); 
            h=symbolparam(nsym,2); 
            
            S1=ceil(w/dx);
            S2=ceil(h/dy); 
            
            if(S1>N1 || S2>N2)
                disp('dimensions of object exceed the dimensions of the system'); 
                return; 
            end
            
            symb=ones(S1,S2); 
            
       end
       if(symboltype(nsym,1)==1)  %this makes circles 
           
            r=symbolparam(nsym,1); 
            
            S1=ceil(2*r/dx);
            S2=ceil(2*r/dy); 
            
            if(S1>N1 || S2>N2)
                disp('dimensions of object exceed the dimensions of the system'); 
                return; 
            end
            
            ra=0:dx:dx*(S1-1); 
            rb=0:dy:dy*(S2-1); 
            
            symb=zeros(S1,S2); 
            %now we need to render the circle
            for na=1:1:S1
                for nb=1:1:S2
                    d=((ra(1,na)-ra(1,end)/2)^2+(rb(1,nb)-rb(1,end)/2)^2)^(1/2); 
                    if(d<r)
                        symb(na,nb)=1; 
                    end
                end
            end
            
            %return; 
            
        end  
    
            
        for nrep=1:1:symbolrep(nsym)
                   
                coord(1,:)=symbolcoord(nsym,:)+symbolvec(nsym,:)*(nrep-1); %we create repeated copies of symbol 
                
                %we need to convert this into absolute pixels for our
                %object, we give the coordinates of one corner
                
                C1=ceil(coord(1,1)/dx);
                C2=ceil(coord(1,2)/dy);
                C1=C1-ceil(S1/2);
                C2=C2-ceil(S2/2);
                
                %first we need to correct C1 and C2 for periodic boundary
                %conditions
                
                C1=C1+1;
                C2=C2+1; 
                
                if(C1>N1)
                    C1=C1-N1; 
                end
                if(C1<=0)
                    C1=C1+N1; 
                end
                if(C2>N2)
                    C2=C2-N2; 
                end
                if(C2<=0)
                    C2=C2+N2; 
                end
                
                %need to calculate mask so we can extract the components of
                %the tensor that are already specified
                
                mask=ones(S1,S2)-symb; 
                master=zeros(S1,S2); 
                
                %specify the contribution of this shape to all of the
                %components of the tensor
                for ntensor=1:1:4
                    for n1=1:1:3
                        for n2=1:1:3
                            master=symbolten(nsym,ntensor,n1,n2)*symb(:,:);
                            if( ((C1+S1-1)<=N1) && ((C2+S2-1)<=N2) )
                                master(:,:)=add(C1:(S1+C1-1),C2:(S2+C2-1),n1,n2,ntensor).*mask(:,:)+master(:,:);
                                add(C1:(S1+C1-1),C2:(S2+C2-1),n1,n2,ntensor)=master(:,:); 
                            end
                            if( ((C1+S1-1)<=N1) && ((C2+S2-1)>N2) )
                                %we need to wrap around part of the tensor
                                master(:,1:(N2-C2+1))=add(C1:(S1+C1-1),C2:N2,n1,n2,ntensor).*mask(:,1:(N2-C2+1))+master(:,1:(N2-C2+1));
                                master(:,(N2-C2+2):S2)=add(C1:(S1+C1-1),1:(S2+C2-N2-1),n1,n2,ntensor).*mask(:,(N2-C2+2):S2)+master(:,(N2-C2+2):S2);
                                add(C1:(S1+C1-1),C2:N2,n1,n2,ntensor)=master(:,1:(N2-C2+1));
                                add(C1:(S1+C1-1),1:(S2+C2-N2-1),n1,n2,ntensor)=master(:,(N2-C2+2):S2);
                            end
                            if( ((C1+S1-1)>N1) && ((C2+S2-1)<=N2) )
                                %we need to wrap around part of the tensor
                                master(1:(N1-C1+1),:)=add(C1:N1,C2:(S2+C2-1),n1,n2,ntensor).*mask(1:(N1-C1+1),:)+master(1:(N1-C1+1),:);
                                master((N1-C1+2):S1,:)=add(1:(S1+C1-N1-1),C2:(S2+C2-1),n1,n2,ntensor).*mask((N1-C1+2):S1,:)+master((N1-C1+2):S1,:);
                                add(C1:N1,C2:(S2+C2-1),n1,n2,ntensor)=master(1:(N1-C1+1),:);
                                add(1:(S1+C1-N1-1),C2:(S2+C2-1),n1,n2,ntensor)=master((N1-C1+2):S1,:);
                            end
                            if( ((C1+S1-1)>N1) && ((C2+S2-1)>N2) )
                                %we need to wrap around part of the tensor
                                master(1:(N1-C1+1),1:(N2-C2+1))=add(C1:N1,C2:N2,n1,n2,ntensor).*mask(1:(N1-C1+1),1:(N2-C2+1))+master(1:(N1-C1+1),1:(N2-C2+1));
                                master(1:(N1-C1+1),(N2-C2+2):S2)=add(C1:N1,1:(S2+C2-N2-1),n1,n2,ntensor).*mask(1:(N1-C1+1),(N2-C2+2):S2)+master(1:(N1-C1+1),(N2-C2+2):S2);
                                master((N1-C1+2):S1,1:(N2-C2+1))=add(1:(S1+C1-N1-1),C2:N2,n1,n2,ntensor).*mask((N1-C1+2):S1,1:(N2-C2+1))+master((N1-C1+2):S1,1:(N2-C2+1));
                                master((N1-C1+2):S1,(N2-C2+2):S2)=add(1:(S1+C1-N1-1),1:(S2+C2-N2-1),n1,n2,ntensor).*mask((N1-C1+2):S1,(N2-C2+2):S2)+master((N1-C1+2):S1,(N2-C2+2):S2);
                                add(C1:N1,C2:N2,n1,n2,ntensor)=master(1:(N1-C1+1),1:(N2-C2+1));
                                add(C1:N1,1:(S2+C2-N2-1),n1,n2,ntensor)=master(1:(N1-C1+1),(N2-C2+2):S2);
                                add(1:(S1+C1-N1-1),C2:N2,n1,n2,ntensor)=master((N1-C1+2):S1,1:(N2-C2+1));
                                add(1:(S1+C1-N1-1),1:(S2+C2-N2-1),n1,n2,ntensor)=master((N1-C1+2):S1,(N2-C2+2):S2);
                            end
                        end
                    end
                end 
                
                %we also want to determine the geometry of the interfaces,
                %all we need to do is add in symbol matrices
                
                if( ((C1+S1-1)<=N1) && ((C2+S2-1)<=N2) )
                    geom(C1:(S1+C1-1),C2:(S2+C2-1))=symb(:,:)+geom(C1:(S1+C1-1),C2:(S2+C2-1)); 
                end
                if( ((C1+S1-1)<=N1) && ((C2+S2-1)>N2) )
                    geom(C1:(S1+C1-1),C2:N2)=symb(:,1:(N2-C2+1))+geom(C1:(S1+C1-1),C2:N2);
                    geom(C1:(S1+C1-1),1:(S2+C2-N2-1))=symb(:,(N2-C2+2):S2)+geom(C1:(S1+C1-1),1:(S2+C2-N2-1));
                end
                if( ((C1+S1-1)>N1) && ((C2+S2-1)<=N2) )
                    geom(C1:N1,C2:(S2+C2-1))=symb(1:(N1-C1+1),:)+geom(C1:N1,C2:(S2+C2-1));
                    geom(1:(S1+C1-N1-1),C2:(S2+C2-1))=symb((N1-C1+2):S1,:)+geom(1:(S1+C1-N1-1),C2:(S2+C2-1));
                end
                if( ((C1+S1-1)>N1) && ((C2+S2-1)>N2) )
                    geom(C1:N1,C2:N2)=symb(1:(N1-C1+1),1:(N2-C2+1))+geom(C1:N1,C2:N2);
                    geom(C1:N1,1:(S2+C2-N2-1))=symb(1:(N1-C1+1),(N2-C2+2):S2)+geom(C1:N1,1:(S2+C2-N2-1));
                    geom(1:(S1+C1-N1-1),C2:N2)=symb((N1-C1+2):S1,1:(N2-C2+1))+geom(1:(S1+C1-N1-1),C2:N2);
                    geom(1:(S1+C1-N1-1),1:(S2+C2-N2-1))=symb((N1-C1+2):S1,(N2-C2+2):S2)+geom(1:(S1+C1-N1-1),1:(S2+C2-N2-1));
                end
                
        end      
 end
 
 toc; 
 
 %for final step we need to extract the tensor components from add
 
 epsilonexp(:,:,:,:)=add(:,:,:,:,1);
 muexp(:,:,:,:)=add(:,:,:,:,2);
 xiexp(:,:,:,:)=add(:,:,:,:,3);
 zetaexp(:,:,:,:)=add(:,:,:,:,4);
 
 
 %we need to provide a set of coordinates which describe the shifts, gives
 %the number of half-pixel shifts
 
 e=zeros(3,3,2);
 e(1,1,:)=[1,0];
 e(2,2,:)=[0,1];
 e(3,3,:)=[0,0];
 e(2,3,:)=[0,0];
 e(3,2,:)=[0,0];
 e(3,1,:)=[0,0];
 e(1,3,:)=[0,0];
 e(1,2,:)=[0,0];
 e(2,1,:)=[0,0];
 
 u=zeros(3,3,2); 
 u(1,1,:)=[0,1];
 u(2,2,:)=[1,0];
 u(3,3,:)=[1,1];
 u(2,3,:)=[1,1];
 u(3,2,:)=[1,1];
 u(3,1,:)=[1,1];
 u(1,3,:)=[1,1];
 u(1,2,:)=[1,1];
 u(2,1,:)=[1,1];
 
 epsilonshift=zeros(N1,N2,3,3,2,2); 
 mushift=zeros(N1,N2,3,3,2,2); 
 
 %need to calculate gradients
 gradx(:,:,1,1)=circshift(geom(:,:),[1,0])-circshift(geom(:,:),[-1,0]); 
 grady(:,:,1,1)=circshift(geom(:,:),[0,1])-circshift(geom(:,:),[0,-1]); 
 
 %we create shifted versions of the geometry and each tensor 
 for n1=0:1:1
     for n2=0:1:1
        if(subpixel~=1) 
            gradx(:,:,n1+1,n2+1)=circshift(gradx(:,:,1,1),[-subpixel/2*n1,-subpixel/2*n2]); 
            grady(:,:,n1+1,n2+1)=circshift(grady(:,:,1,1),[-subpixel/2*n1,-subpixel/2*n2]);
            epsilonshift(:,:,:,:,n1+1,n2+1)=circshift(epsilonexp(:,:,:,:),[-subpixel/2*n1,-subpixel/2*n2]); 
            mushift(:,:,:,:,n1+1,n2+1)=circshift(muexp(:,:,:,:),[-subpixel/2*n1,-subpixel/2*n2]);
        else
            gradx(:,:,n1+1,n2+1)=circshift(gradx(:,:,1,1),[0/2*n1,0/2*n2]); 
            grady(:,:,n1+1,n2+1)=circshift(grady(:,:,1,1),[0/2*n1,0/2*n2]);
            epsilonshift(:,:,:,:,n1+1,n2+1)=circshift(epsilonexp(:,:,:,:),[0/2*n1,0/2*n2]); 
            mushift(:,:,:,:,n1+1,n2+1)=circshift(muexp(:,:,:,:),[0/2*n1,0/2*n2]);
        end
     end
 end
 
 
%now we need to compress dimensions
 
N1=sx*Repeatx;
N2=sy*Repeaty; 
N3=sz*Repeatz; 

Nxyz=N1*N2*N3; 
 
epsilon=zeros(N1,N2,3,3); 
mu=zeros(N1,N2,3,3);
xi=zeros(N1,N2,3,3); 
zeta=zeros(N1,N2,3,3);  

avggradx=zeros(N1,N2,2,2); 
avggrady=zeros(N1,N2,2,2); 
edge=sparse(N1,N2); 

tic; 


%we compute average of gradient in each of the 4 geometry grids we need to
%maintain

for na=1:1:2
    for nb=1:1:2
        for n1=1:1:N1
            for n2=1:1:N2   
                avggradx(n1,n2,na,nb)=mean(mean(gradx(((n1-1)*subpixel+1):(n1*subpixel),((n2-1)*subpixel+1):(n2*subpixel),na,nb),1),2); 
                avggrady(n1,n2,na,nb)=mean(mean(grady(((n1-1)*subpixel+1):(n1*subpixel),((n2-1)*subpixel+1):(n2*subpixel),na,nb),1),2);      
            end
        end
    end
end

%we just downsample xi and zeta, we don't worry about the grids
for n1=1:1:N1
    for n2=1:1:N2   
        xi(n1,n2,:,:)=xiexp((n1*subpixel),(n2*subpixel),:,:);
        zeta(n1,n2,:,:)=zetaexp((n1*subpixel),(n2*subpixel),:,:);
    end
end

 %we downsample everything, depending on the tensor component, we
 %downsample using a different grid, specify which grid with e and u
 for na=1:1:3
     for nb=1:1:3
         for n1=1:1:N1
             for n2=1:1:N2
                epsilon(n1,n2,na,nb)=epsilonshift((n1*subpixel),(n2*subpixel),na,nb,e(na,nb,1)+1,e(na,nb,2)+1);
                mu(n1,n2,na,nb)=mushift((n1*subpixel),(n2*subpixel),na,nb,u(na,nb,1)+1,u(na,nb,2)+1);  
             end
         end
     end
 end
        


toc; 

tic; 

%need to do seperate loops for mu and epsilon

epsilon=subpixelaverage(avggradx,avggrady,e,epsilonshift,epsilon,subpixel);

mu=subpixelaverage(avggradx,avggrady,u,mushift,mu,subpixel);



toc; 

%%%%%%%%display critical tensor components

%{

figure; 
imagesc(r2(:,2),r1(:,1),reshape(epsilon(:,:,1,1),N1,N2)); 
axis equal; 

figure; 
imagesc(r2(:,2),r1(:,1),reshape(epsilon(:,:,3,3),N1,N2)); 
axis equal; 

figure; 
imagesc(r2(:,2),r1(:,1),reshape(mu(:,:,1,1),N1,N2)); 
axis equal; 

figure; 
imagesc(r2(:,2),r1(:,1),abs(reshape(xi(:,:,1,2),N1,N2))); 
axis equal; 

figure; 
imagesc(r2(:,2),r1(:,1),abs(reshape(xi(:,:,1,1),N1,N2))); 
axis equal; 

return; 
%}


%only the first splitnum file outputs the system information, redundant to
%do otherwise
if(splitnum==1)

    %we need to give sufficient information so that we can open all files
    %automatically

    fp=fopen([dirt,'/pos.txt'],'w');  %print out position information
    fprintf(fp,'%d \n',N1);
    fprintf(fp,'%d \n',N2); 
    for n=1:1:N1
        fprintf(fp,'%d \n',norm(r1(n,:))); 
    end
    for n=1:1:N2
        fprintf(fp,'%d \n',norm(r2(n,:))); 
    end
    fclose(fp); 
    
    fp=fopen([dirt,'/epsilonreal.txt'],'w'); 
    for n1=1:1:N1
        for n2=1:1:N2
            for na=1:1:3
                for nb=1:1:3
                    fprintf(fp,'%d ',real(epsilon(n1,n2,na,nb))); 
                end
            end
            fprintf(fp,' \n'); 
        end
    end
    fclose(fp); 
 
    fp=fopen([dirt,'/epsilonimag.txt'],'w'); 
    for n1=1:1:N1
        for n2=1:1:N2
            for na=1:1:3
                for nb=1:1:3
                    fprintf(fp,'%d ',imag(epsilon(n1,n2,na,nb))); 
                end
            end
            fprintf(fp,' \n'); 
        end
    end
    fclose(fp); 
    
    fp=fopen([dirt,'/mureal.txt'],'w'); 
    for n1=1:1:N1
        for n2=1:1:N2
            for na=1:1:3
                for nb=1:1:3
                    fprintf(fp,'%d ',real(mu(n1,n2,na,nb))); 
                end
            end
            fprintf(fp,' \n'); 
        end
    end
    fclose(fp); 
    
     fp=fopen([dirt,'/muimag.txt'],'w'); 
    for n1=1:1:N1
        for n2=1:1:N2
            for na=1:1:3
                for nb=1:1:3
                    fprintf(fp,'%d ',imag(mu(n1,n2,na,nb))); 
                end
            end
            fprintf(fp,' \n'); 
        end
    end
    fclose(fp); 
    
    fp=fopen([dirt,'/xireal.txt'],'w'); 
    for n1=1:1:N1
        for n2=1:1:N2
            for na=1:1:3
                for nb=1:1:3
                    fprintf(fp,'%d ',real(xi(n1,n2,na,nb))); 
                end
            end
            fprintf(fp,' \n'); 
        end
    end
    fclose(fp); 
    
    fp=fopen([dirt,'/xiimag.txt'],'w'); 
    for n1=1:1:N1
        for n2=1:1:N2
            for na=1:1:3
                for nb=1:1:3
                    fprintf(fp,'%d ',imag(xi(n1,n2,na,nb))); 
                end
            end
            fprintf(fp,' \n'); 
        end
    end
    fclose(fp); 
    
    fp=fopen([dirt,'/zetareal.txt'],'w'); 
    for n1=1:1:N1
        for n2=1:1:N2
            for na=1:1:3
                for nb=1:1:3
                    fprintf(fp,'%d ',real(zeta(n1,n2,na,nb))); 
                end
            end
            fprintf(fp,' \n'); 
        end
    end
    fclose(fp); 
 
    fp=fopen([dirt,'/zetaimag.txt'],'w'); 
    for n1=1:1:N1
        for n2=1:1:N2
            for na=1:1:3
                for nb=1:1:3
                    fprintf(fp,'%d ',imag(zeta(n1,n2,na,nb))); 
                end
            end
            fprintf(fp,' \n'); 
        end
    end
    fclose(fp); 
 
    

end



%we need to reform our tensors for the script written to solve the system
%in real space

epsilonp=zeros(1,N1*N2*N3,3,3); 
mup=zeros(1,N1*N2*N3,3,3); 
xip=zeros(1,N1*N2*N3,3,3); 
zetap=zeros(1,N1*N2*N3,3,3); 
fill=zeros(1,N1*N2*N3); 

%we need to make sure we are reforming the rows in the correct way

for n1=1:1:3
    for n2=1:1:3
        count=0; 
        for nb=1:1:N2
            for na=1:1:N1
                count=count+1;
                epsp(1,count,n1,n2)=epsilon(na,nb,n1,n2); 
                mup(1,count,n1,n2)=mu(na,nb,n1,n2); 
                xip(1,count,n1,n2)=xi(na,nb,n1,n2); 
                zetap(1,count,n1,n2)=zeta(na,nb,n1,n2);
            end
        end
    end
end


%% 

%disp(epsp(1,:,1,1)); 
%return; 

%%%now we need to put into B1 matrix to be used in external code

B=spdiags([zetap(1,:,3,1),fill,fill,fill,fill,fill;
    zetap(1,:,2,1),zetap(1,:,3,2),fill,fill,fill,fill;
    zetap(1,:,1,1),zetap(1,:,2,2),zetap(1,:,3,3),fill,fill,fill;
    epsp(1,:,3,1),zetap(1,:,1,2),zetap(1,:,2,3),mup(1,:,3,1),fill,fill;
    epsp(1,:,2,1),epsp(1,:,3,2),zetap(1,:,1,3),mup(1,:,2,1),mup(1,:,3,2),fill;
    epsp(1,:,1,1),epsp(1,:,2,2),epsp(1,:,3,3),mup(1,:,1,1),mup(1,:,2,2),mup(1,:,3,3);
    fill,epsp(1,:,1,2),epsp(1,:,2,3),xip(1,:,3,1),mup(1,:,1,2),mup(1,:,2,3);
    fill,fill,epsp(1,:,1,3),xip(1,:,2,1),xip(1,:,3,2),mup(1,:,1,3);
    fill,fill,fill,xip(1,:,1,1),xip(1,:,2,2),xip(1,:,3,3);
    fill,fill,fill,fill,xip(1,:,1,2),xip(1,:,2,3);
    fill,fill,fill,fill,fill,xip(1,:,1,3)].',[-5*Nxyz,-4*Nxyz,-3*Nxyz,-2*Nxyz,-Nxyz,0,Nxyz,2*Nxyz,3*Nxyz,4*Nxyz,5*Nxyz],6*Nxyz,6*Nxyz);

%we want to reduce to a subblock so we can solve for TE/TM bandstructure
B=circshift(B,[N1*N2*N3,N1*N2*N3]); 
B=B((3*N1*N2*N3+1):(6*N1*N2*N3),(3*N1*N2*N3+1):(6*N1*N2*N3));

%if we calculate inner products we need the B matrix
if(splitnum==1)
   
    [ids,jds,val]=find(B); 
    
    si=size(ids);
    si=si(1,1); 
    
    fp=fopen([dirt,'/Breal.txt'],'w'); 
    fprintf(fp,'%d \n',si); 
    for n1=1:1:si
        fprintf(fp,'%d %d %d \n',ids(n1,1),jds(n1,1),real(val(n1,1))); 
    end
    fclose(fp); 
    
    fp=fopen([dirt,'/Bimag.txt'],'w'); 
    fprintf(fp,'%d \n',si); 
    for n1=1:1:si
        fprintf(fp,'%d %d %d \n',ids(n1,1),jds(n1,1),imag(val(n1,1))); 
    end
    fclose(fp); 
    
end


%matlabpool open;
%matlabpool(4);

s=size(k); 
sk=s(1,2); 

data=zeros(sk,bands); 
dx=norm(R1)*Repeatx/N1; 
dy=norm(R2)*Repeaty/N2; 
dz=norm(R1)*Repeatx/N1; %essentially garabage value we put in as a placeholder

opts.tol=10*10^-14;
opts.maxit=30000; 


%first circshift B into the appropriate position




for ind = (1+ceil(sk*(splitnum-1)/split)):1:ceil(sk*splitnum/split)
    kpoint=k(1,ind)*G1+k(2,ind)*G2; 
    kx = kpoint(1,1); 
    ky = kpoint(1,2); 
    kz = kpoint(1,3); 
    disp('creating matrix');
    A=curlperiodicv2new(N1,N2,N3,dx,dy,dz,kx,ky,kz,Mx,My,Mz);
    
    A=circshift(A,[N1*N2*N3,N1*N2*N3]); 
    A=A((3*N1*N2*N3+1):(6*N1*N2*N3),(3*N1*N2*N3+1):(6*N1*N2*N3));
    
    disp('running eigs');
    tic; 
    
    [V,D]=eigs(A,B,bands,center,opts); 
    
    toc; 
    disp('finished eigs'); 
    
    %extract eigenvalues from diagonal matrix elements
    omg=diag(D)/(2*pi); 
    omgs=sort(omg);
    
    disp('writing data');
    disp(ind); 
    a=sprintf('%d',ind); 
    fname=[dirt,'/','output',a,'.txt']; 
    fp=fopen(fname,'w');
    for n=1:1:bands
        fprintf(fp,'%d \n',omgs(n,1));
    end
    fclose(fp);
    
    sd=size(kind); 
    sd=sd(1,2); 
    
    %we want to sort eigenvalues and eigenvectors before outputting will
    %salve us a lot of time 
   
    locsort=cell(bands,2);
    
    for nd=1:1:bands
        locsort(nd,1)={real(omg(nd,1))};
        locsort(nd,2)={V(:,nd)}; 
    end
        
    disp(omg(:,1)); 
    disp(locsort); 
    
    sorted=sortrows(locsort,1); %sort based on omega values
        
    
    
    for nd=1:1:bands
        omg(nd,1)=cell2mat(sorted(nd,1)); 
        V(:,nd)=cell2mat(sorted(nd,2));  
    end
    
    
    for nd=1:1:sd
    
        if(ind==kind(1,nd))  %we only want to pick out the surface states at a particular k-point
    
           
            fname=[dirt,'/','eigenvalue',a,'.txt']; 
            fp2=fopen(fname,'w'); 
        
            %we iterate through the bands we found and we see if any of them fall
            %in the range where we want to keep the surface states
            veccount=1;
            for n=1:1:bands
                if(omg(n,1)<omegamax && omg(n,1)>omegamin)
                    disp('output eigenvector'); 
                    fprintf(fp2,'%d \n',omg(n,1));
                    disp(omg(n,1)); 
                    %disp(V(:,n));
                    b=sprintf('%d',veccount); 
                    
                    fname=[dirt,'/','eigenvectorreal',a,'vec',b,'.txt']; 
                    fp0=fopen(fname,'w'); 
                    fname=[dirt,'/','eigenvectorimag',a,'vec',b,'.txt']; 
                    fp1=fopen(fname,'w'); 
                    
                    for n2=1:1:(N1*N2*N3*3)
                        fprintf(fp0,'%d ',real(V(n2,n))); 
                        fprintf(fp1,'%d ',imag(V(n2,n))); 
                    end
                    
                    fclose(fp0); 
                    fclose(fp1);
                    
                    veccount=veccount+1;
                end
            end
        
            fclose(fp2); 
        end
    
    end
    
    
    
    
end



return; 


