function [gamma] = kaewu_etran(coordi, coordj, webdir,L)
% returns the rotation matrix gamma given input L from other function%


fr = [((coordj(1)-coordi(1))/L), ((coordj(2)-coordi(2))/L), ((coordj(3)-coordi(3))/L)]; 
sr = [webdir(1),webdir(2),webdir(3)];
tr = (cross(fr,sr)); 


gamma = [ fr 0 0 0 0 0 0 0 0 0;
          sr 0 0 0 0 0 0 0 0 0;
          tr 0 0 0 0 0 0 0 0 0;
          0 0 0 fr 0 0 0 0 0 0;
          0 0 0 sr 0 0 0 0 0 0;  
          0 0 0 tr 0 0 0 0 0 0;
          0 0 0 0 0 0 fr 0 0 0; 
          0 0 0 0 0 0 sr 0 0 0;
          0 0 0 0 0 0 tr 0 0 0;
          0 0 0 0 0 0 0 0 0 fr;
          0 0 0 0 0 0 0 0 0 sr;
          0 0 0 0 0 0 0 0 0 tr;];
     

end

