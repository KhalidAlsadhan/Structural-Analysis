function [FeF] = kaewu_computeMemberFEFs(w,L)

% compute FeF using the same DOF numbering as memb_id
%   dx1, dy1, dz1, thetax1, thetay1, thetaz1, ...
%   dx2, dy2, dz2, thetax2, thetay2, thetaz2

% print out distributed load to check
w ; 

% x-component distributed load
w_x = w(1);
% y-component distributed load
w_y = w(2);
% z-component distributed load
w_z = w(3);


% Rx due to w_x; Mx is 0 since x along element
Rx = (w_x*L/2);
Mx = 0;

% Ry due to w_y; My due to w_z
Ry = (w_y*L/2);
My = w_z*L^2/12;

% Rz due to w_z; Mz due to w_y
Rz = (w_z*L/2);
Mz = w_y*L^2/12;

% node i has reaction moments opposite in sign of load 
% node j has reaction moments same sign of load if 
% counterclockwise positive and tensile loads positive
FeF_i = [-Rx, -Ry, -Rz, 0, -My, -Mz];
FeF_j = [-Rx, -Ry, -Rz, 0,  My,  Mz];

FeF = [ FeF_i, FeF_j];


% note that FeF have a positive rotation at the 2nd node of a member

   
end

