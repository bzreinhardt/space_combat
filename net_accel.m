function [net_accel] = net_accel(x,y,z,x1,thrust_ratio,A,solar_power,m)
%Net accel finds the magnitude and sign of the net acceleration caused
%by gravity and station keeping thrusters.  If the acceleration is positive
%the thrusters cannot compensate for the gravitational acceleration
r_from_sun = ((x-x1)^2+y^2+z^2)^.5;
L3 = fzero(@x_accel,-1E11);
flux = solar_power/(4*pi*r_from_sun^2);

thrust = thrust_ratio*flux*A;
thruster_accel = thrust/m; 

grav_accel = norm(accel(x,y,z),2);

net_accel = grav_accel-thruster_accel; 



end

