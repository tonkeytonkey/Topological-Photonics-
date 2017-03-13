
close all; 
clear all; 
clc; 


sx=32;
sy=32; 
sz=1; 

subpixel=1; 


dirt=['rodsweep1/sweep1']; 

I=i; 

fp=fopen([dirt,'/settings1.txt'],'r');  %print out position information
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
  
      if(symboltype(nsym,1)==2 )  %this makes blocks 
           
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


%we want to construct shifted versions of geom so we can calculate the area
%exactly 




%we want to compute the mean of the different geometry grids so we can get
%the area of each tensor in each region

%{

figure; 
imagesc(r2(:,2),r1(:,1),reshape(avggeom(:,:,1,1),N1,N2)); 
axis equal; 
title('avggeom 1 1 ');

figure; 
imagesc(r2(:,2),r1(:,1),reshape(avggeom(:,:,2,1),N1,N2)); 
axis equal; 
title('avggeom 2 1 ');

figure; 
imagesc(r2(:,2),r1(:,1),reshape(avggeom(:,:,1,2),N1,N2)); 
axis equal; 
title('avggeom 1 2 ');

%}



for na=1:1:2
    for nb=1:1:2
    
        %{
figure; 
imagesc(r2(:,2),r1(:,1),reshape(gradx(:,:,na,nb),N1*subpixel,N2*subpixel)); 
axis equal; 
title('gradx');
%}
        
%{
figure; 
imagesc(r2(:,2),r1(:,1),reshape(grady(:,:,na,nb),N1*subpixel,N2*subpixel)); 
axis equal; 
title('grady');
%}
       
        %{
figure; 
imagesc(r2(:,2),r1(:,1),reshape(avggradx(:,:,na,nb),N1,N2)); 
axis equal; 
title('avggradx');
%}
        
%{
figure; 
imagesc(r2(:,2),r1(:,1),reshape(avggrady(:,:,na,nb),N1,N2)); 
axis equal; 
title('avggrady');
%}

    end
end


%we want to mask out the sections associated with the gradient and show
%that it just contains the pixels we averaged

mask=zeros(N1,N2); 
for n1=1:1:N1
    for n2=1:1:N2
        if( (abs(avggradx(n1,n2,2,1))+abs((avggrady(n1,n2,2,1))))>0 )
            mask(n1,n2)=epsilon(n1,n2,1,1); 
        end
    end
end




%now calculate material tensors

%{
figure; 
imagesc(r2(:,2),r1(:,1),reshape(epsilon(:,:,1,1),N1,N2)); 
axis equal; 
title('Epsilon 1,1 real');

figure; 
imagesc(r2(:,2),r1(:,1),reshape(epsilon(:,:,1,2),N1,N2)); 
axis equal; 
title('Epsilon 1,2 real');

figure; 
imagesc(r2(:,2),r1(:,1),reshape(epsilon(:,:,2,2),N1,N2)); 
axis equal; 
title('Epsilon 2,2 real');

figure; 
imagesc(r2(:,2),r1(:,1),reshape(epsilon(:,:,3,1),N1,N2)); 
axis equal; 
title('Epsilon 3,1 real');
%}

figure; 
imagesc(r2(:,2),r1(:,1),reshape(epsilon(:,:,1,1),N1,N2)); 
axis equal; 
title('epsilon 1,1 real'); 

figure; 
imagesc(r2(:,2),r1(:,1),reshape(imag(mu(:,:,1,2)),N1,N2)); 
axis equal; 
title('mu 1,2 imag'); 

figure; 
imagesc(r2(:,2),r1(:,1),reshape(epsilon(:,:,3,3),N1,N2)); 
axis equal; 
title('epsilon 3,3 real'); 


%{
figure; 
imagesc(r2(:,2),r1(:,1),reshape(epsilon(:,:,3,1),N1,N2)); 
axis equal; 
title('epsilon 1,3 real'); 
%}

%{
figure; 
imagesc(r2(:,2),r1(:,1),reshape(mask(:,:),N1,N2)); 
axis equal; 
title('epsilon 3,3 real with overlay'); 
%}


figure; 
imagesc(r2(:,2),r1(:,1),reshape(mu(:,:,3,3),N1,N2)); 
axis equal; 
title('mu 3,3 real'); 


%{

figure; 
imagesc(r2(:,2),r1(:,1),reshape(mu(:,:,2,2),N1,N2)); 
axis equal; 
title('mu 2,2 real'); 

figure; 
imagesc(r2(:,2),r1(:,1),reshape(mu(:,:,3,3),N1,N2)); 
axis equal; 
title('mu 3,3 real'); 



figure; 
imagesc(r2(:,2),r1(:,1),real(reshape(xi(:,:,1,3),N1,N2))); 
axis equal; 
title('xi 3,1 real'); 

figure; 
imagesc(r2(:,2),r1(:,1),real(reshape(xi(:,:,3,1),N1,N2))); 
axis equal; 
title('xi 1,3 real'); 
%}

%{
figure; 
imagesc(r2(:,2),r1(:,1),real(reshape(xi(:,:,1,1),N1,N2))); 
axis equal; 
title('xi 1 1 real'); 
%}



return; 
