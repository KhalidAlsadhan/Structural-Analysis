function [elk] = kaewu_estiff(A, Izz, Iyy, J, Ayy, Azz, E, v, L)

G = (E/(2*(1+v))); 
n = ((Izz)/(A*(5/6)*G));
coff = ((Izz)/(L^3/(12) +L*n));
%The Above variables are used to calculate the shear coefficient 
%to change the affected indices. 

   
   
elk =            E*[(A/L) 0 0 0 0 0 -(A/L) 0 0 0 0 0;
                  0 coff 0 0 0 (coff*L/2) 0 -(coff) 0 0 0 (coff*L/2);
                  0 0 (12*Iyy/L^3) 0 -(6*Iyy*1/L^2) 0 0 0 -(12*Iyy/L^3) 0 -(6*Iyy/L^2) 0;
                  0 0 0 (J/(2*(1+v)*L)) 0 0 0 0 0 -(J/(2*(1+v)*L)) 0 0 ;
                  0 0 -(6*Iyy*1/L^2) 0 (4*Iyy/L) 0 0 0 -(6*Iyy*1/L^2) 0 (2*Iyy/L) 0 ;
                  0 (coff*L/2) 0 0 0 (coff*((L^2)/3 +n)) 0 -((coff*L/2)) 0 0 0 coff*((L^2)/6 - n);
                  -(A/L) 0 0 0 0 0 (A/L) 0 0 0 0 0;
                  0 -(coff) 0 0 0 -(coff*L/2) 0 (coff) 0 0 0 -(coff*L/2);
                  0 0 -(12*Iyy/L^3) 0 (6*Iyy*1/L^2) 0 0 0 (12*Iyy/L^3) 0 (6*Iyy*1/L^2) 0;
                  0 0 0 -(J/(2*(1+v)*L)) 0 0 0 0 0 (J/(2*(1+v)*L)) 0 0 ;
                  0 0 -(6*Iyy*1/L^2) 0 (2*Iyy/L) 0 0 0 (6*Iyy*1/L^2) 0 (4*Iyy/L) 0; 
                  0 (coff*L/2) 0 0 0 coff*(L^2/6 - n) 0 -((coff*L/2)) 0 0 0 coff*(L^2/3 +n)] ;
              
              
       


end