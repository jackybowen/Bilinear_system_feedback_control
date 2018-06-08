function [ dydt ] = f_1d_burgers( t, y )
% The 1-D Burgers' equation by the Finite Difference Method (a time march)
% Numerical scheme used is a first order upwind in time
...and space for the convection terms 
...and a second order central difference in space
...for the diffusive term
global dx
x = 0:dx:2;
y
for i = 1:length(x)
    y(i)
end


end

