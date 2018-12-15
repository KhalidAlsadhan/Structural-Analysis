% length function
function [L] = lengthfunction(coordi, coordj)

%this function calculates the length given the coordinates of node i and j
L = sqrt(((coordi(1)-coordj(1))^2 + (coordi(2)-coordj(2))^2 +(coordj(3)-coordi(3))^2));


end