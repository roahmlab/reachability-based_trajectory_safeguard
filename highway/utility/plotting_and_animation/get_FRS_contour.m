function [FRS_contour]=get_FRS_contour(FRS,w_l,z0,k_in)
           
           
            
          
                 Dx=FRS.Dx;
                 Dy=FRS.Dy;
                k=FRS.k;
                z=FRS.z;
                wk=msubs(w_l,k,k_in);
                 x0=FRS.x0;
                 y0=FRS.y0;
            %plot(z0(1),z0(2),'*')
             xvec = linspace(-1,1,100) ;
             yvec = linspace(-1,1,100) ;
             [X,Y] = meshgrid(xvec, yvec) ;
             ZZ = [X(:) Y(:)]' ;
             F = reshape(full(msubs(wk,z(1:2), ZZ)),100,100) ;

             [hh]=contourc(xvec,yvec,F,[0 0]);
              
               idxs=find(hh(2,:)>=2);
               [idm,I]=max(idxs);
               if length(idxs)-I>=1
                   x=hh(1,idm+1:idxs(I+1)-1);
                   y=hh(1,idm+1:idxs(I+1)-1);
               else
                   x=hh(1,idm+1:end);
                   y=hh(2,idm+1:end);
               end
               
                  xshift=(x-x0/Dx)*cos(z0(3))-sin(z0(3))*(y-y0/Dy);
                  yshift=(y-y0/Dy)*cos(z0(3))+sin(z0(3))*(x-x0/Dx);
               
               
              FRS_contour=z0(1:2)+[xshift*Dx;yshift*Dy];
             
            end
           
                
  
  