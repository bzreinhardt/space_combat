function [accel] = tot_accel(x,y,z)
%bunch of constants
G = 6.67300E-11;
m_E = 5.9742E24; %kg
m_S = 1.98892E30; %kg
r_E = 149598000E3;%meters
L3 = fzero(@x_accel,-1E11);

mu = m_E/(m_E+m_S);

M1 = m_S;
M2 = mu*M1/(1-mu);



a = r_E;
x2 = a;
x1 = -x2*M2/M1;
r1 = ((x-x1)^2+y^2+z^2)^.5;
r2 = ((x-x2).^2+y^2+z^2).^.5;


n = (G*(M1+M2)/x2^3)^.5;
accel_y = n^2*y - G*M1/(r1^3)*y - G*M2/(r2^3)*y;
accel_x = n^2*x - G*M1/(r1^3)*(x-x1) - G*M2/(r2^3)*(x-x2);
accel_z = G*M1/(r1^3)*z - G*M2/(r2^3)*z;

accel = [accel_x;accel_y;accel_z]; 

end

