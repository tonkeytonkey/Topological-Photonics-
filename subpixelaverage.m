
%we take in tensor epsilon (or mu) and conduct the correct average at the
%surface based on the tensor component

function [epsilon] = subpixelaverage(avggradx,avggrady,e,epsilonshift,epsilon,subpixel)

for na=1:1:3 
    for nb=1:1:3
    

        edge(:,:)=abs(avggradx(:,:,e(na,nb,1)+1,e(na,nb,2)+1))+abs(avggrady(:,:,e(na,nb,1)+1,e(na,nb,2)+1));

        [iedge,jedge,valedge]=find(sparse(edge));

        %now we want to iterate through points where the average gradient is not
        %equal to zero and conduct an effective average of the tensors

        iedge=iedge.'; %transpose so dimensions work out the way we want
        jedge=jedge.'; 

        s=size(iedge);
        sedges=s(1,2); 

        epsloc=zeros(subpixel,subpixel,3,3);
        epseff=zeros(subpixel,subpixel,3,3);
        epsavg=zeros(3,3);
        
        disp(sedges); 

        %iterate through all edges of the system
        for n=1:1:sedges
    
            %first compute rotation matrix
            gx=avggradx(iedge(1,n),jedge(1,n),e(na,nb,1)+1,e(na,nb,2)+1);
            gy=avggrady(iedge(1,n),jedge(1,n),e(na,nb,1)+1,e(na,nb,2)+1);
            ang=acos(gx/(gx^2+gy^2)^(1/2)); 
            R=[cos(ang),-sin(ang),0;sin(ang),cos(ang),0;0,0,1]; 
    
            %compute local epsilon and mu tensors in the cell 
            epsloc=epsilonshift( ((iedge(1,n)-1)*subpixel+1):(iedge(1,n)*subpixel),((jedge(1,n)-1)*subpixel+1):(jedge(1,n)*subpixel),:,:,e(na,nb,1)+1,e(na,nb,2)+1);
    
            %now we need to rotate the tensor for each individual pixel
            for n1=1:1:subpixel
                for n2=1:1:subpixel
                    epsloc(n1,n2,:,:)=R\reshape(epsloc(n1,n2,:,:),3,3)*R;
                end
            end
    
            %now that we've transformed each local point, we need to rearrange the
            %tenors into effective form, average them, and then invert the average
    
            epseff(:,:,1,1)=-1./epsloc(:,:,1,1); 
            epseff(:,:,1,2)=epsloc(:,:,1,2)./epsloc(:,:,1,1); 
            epseff(:,:,1,3)=epsloc(:,:,1,3)./epsloc(:,:,1,1); 
            epseff(:,:,2,1)=epsloc(:,:,2,1)./epsloc(:,:,1,1); 
            epseff(:,:,3,1)=epsloc(:,:,3,1)./epsloc(:,:,1,1); 
            epseff(:,:,2,2)=epsloc(:,:,2,2)-(epsloc(:,:,2,1).*epsloc(:,:,1,2))./epsloc(:,:,1,1); 
            epseff(:,:,3,3)=epsloc(:,:,3,3)-(epsloc(:,:,3,1).*epsloc(:,:,1,3))./epsloc(:,:,1,1); 
            epseff(:,:,2,3)=epsloc(:,:,2,3)-(epsloc(:,:,2,1).*epsloc(:,:,1,3))./epsloc(:,:,1,1); 
            epseff(:,:,3,2)=epsloc(:,:,3,2)-(epsloc(:,:,3,1).*epsloc(:,:,1,2))./epsloc(:,:,1,1); 
    
            %now we need to average this within each cell
    
            epseffavg=mean(mean(epseff(1:subpixel,1:subpixel,:,:),1),2); 
    
            epseffavg=reshape(epseffavg(1,1,:,:),3,3);
    
            %now we need to invert these effective matrices, and then rotate the
            %result back to the proper coordinate system
    
            epsavg(1,1)=-1./epseffavg(1,1); 
            epsavg(1,2)=-epseffavg(1,2)./epseffavg(1,1); 
            epsavg(1,3)=-epseffavg(1,3)./epseffavg(1,1); 
            epsavg(2,1)=-epseffavg(2,1)./epseffavg(1,1); 
            epsavg(3,1)=-epseffavg(3,1)./epseffavg(1,1); 
            epsavg(2,2)=epseffavg(2,2)-(epseffavg(2,1).*epseffavg(1,2))./epseffavg(1,1); 
            epsavg(3,3)=epseffavg(3,3)-(epseffavg(3,1).*epseffavg(1,3))./epseffavg(1,1); 
            epsavg(2,3)=epseffavg(2,3)-(epseffavg(2,1).*epseffavg(1,3))./epseffavg(1,1); 
            epsavg(3,2)=epseffavg(3,2)-(epseffavg(3,1).*epseffavg(1,2))./epseffavg(1,1); 
    
    
            %now we need to unrotate these parts

            epsavg=R*epsavg/R;
    
            
            %finally set the appropriate cell on the tensors
            epsilon(iedge(1,n),jedge(1,n),na,nb)=epsavg(na,nb);  %only pick up the appropriate component
            
        end
    end 
end

return; 