G = 6.67300E-11;
m_E = 5.9742E24; %kg
m_S = 1.98892E30; %kg
r_E = 149598000E3;%meters
R_sun = 6.955E8;

mu = m_E/(m_E+m_S);

M1 = m_S;
M2 = mu*M1/(1-mu);

pwr2thrust = 12/500E3; %N/watt



%% Finding L3
%general measurements
a = r_E;
x2 = a;
x1 = -x2*M2/M1;


n = (G*(M1+M2)/x2^3)^.5;
x = -1*logspace(6,15,10000);
y = 0;
z = 0;
r1 = ((x-x1).^2+y^2+z^2).^.5;
r2 = ((x-x2).^2+y^2+z^2).^.5;

a_x = n^2*x - G*M1./(r1.^3).*(x-x1) - G*M2./(r2.^3).*(x-x2);

figure(1)
semilogx(x,a_x,x,0)
ylim([-10,10])
title('X gravitational acceleration to left of M1')
xlabel('x(m)')
ylabel('x''(m/s^2)')
%fzero
%Accelerations at graviational point
%a_x = n^2*x - G*M1/(r1^3)*(x-x1) - G*M2/(r2^3)*(x-x2);
%a_y = n^2*y - G*M1/(r1^3)*y - G*M2/(r2^3)*y;
%a_z = - G*M1/(r1^3)*z - G*M2/(r2^3)*z;
L3 = fzero(@x_accel,-1E11);
max_y_dev = (-L3+x2)*R_sun/(x2-x1);

perturb_x = linspace(-100000,100000,700);
perturb_y = linspace(0,max_y_dev,700);


circ_y = linspace(0,abs(max_y_dev),5000);
circ_x = -1*((L3-x1)^2-circ_y.^2).^.5+x1;

perturb_x_accel = zeros(length(perturb_x));
perturb_y_accel = zeros(length(perturb_y));
for i = 1:length(perturb_x)
perturb_x_accel(i) = norm(accel(perturb_x(i)+L3,0,0),2);
perturb_y_accel(i) = norm(accel(L3,perturb_y(i),0),2);
end


perturb_circ_accel = zeros(length(circ_y));
angle = zeros(length(circ_y));
for i = 1:length(circ_y)
    perturb_circ_accel(i) = norm(accel(circ_x(i),circ_y(i),0),2);
    angle(i) = atan(circ_y(i)/circ_x(i));
end
figure(2)
clf
subplot(311)
plot(perturb_x,perturb_x_accel);
title('Accelerations resulting from x perturbations from the L3 point')
subplot(312)
plot(perturb_y,perturb_y_accel);
xlabel('perturbation in y direction from L3 (m)');ylabel('acceleration magnitude(m/s^2');
title('Accelerations resulting from y perturbations from the L3 point')
subplot(313)
semilogy(angle*180/pi, perturb_circ_accel);
xlabel('Perturbations along the circle of radius(L3-x1)');
title('Accelerations resulting from circular perturbations from the L3 point')


%%
%Find flux at orbit location
l = 100; w = 20; h = 20;
A = h*l;
rho1 = 280;
rho2 = 2800;
m1 = l*w*h*rho1;
m2 = l*w*h*rho2;

F_E = 1370; %w/m
F_r = F_E * r_E^2/(L3-x1)^2;
power = F_r*A;


thrust_ratio = 5/200000; %thrust ratio of a vasmir thruster in N/Watt
thrust_ratio2 = .5/200000;
solar_power = 3.85285E26; %watts

%create grid around L3
position_x = linspace(L3-10000000,L3+10000000,100);
position_y = linspace(-max_y_dev,max_y_dev,100);
%accel_mag1 = zeros(length(position_y),length(position_x));
%accel_test1 = zeros(length(position_y),length(position_x));
%accel_mag2= zeros(length(position_y),length(position_x));
%accel_test2 = zeros(length(position_y),length(position_x));
accel_mag3= zeros(length(position_y),length(position_x));
accel_test3 = zeros(length(position_y),length(position_x));
for i = 1:length(position_x)
    for j = 1:length(position_y)
       %accel_test1(j,i) = norm(accel(position_x(i),position_y(j),0),2);
        %accel_mag1(j,i) = net_accel(position_x(i),position_y(j),0,x1,thrust_ratio,A,solar_power,m1);
        %accel_test2(j,i) = norm(accel(position_x(i),position_y(j),0),2);
        %accel_mag2(j,i) = net_accel(position_x(i),position_y(j),0,x1,thrust_ratio,A,solar_power,m2);
        accel_test3(j,i) = norm(accel(position_x(i),position_y(j),0),2);
        accel_mag3(j,i) = net_accel(position_x(i),position_y(j),0,x1,thrust_ratio2,A,solar_power,m2);
    end
end

%figure(5)
%surfc(position_x-L3,position_y,accel_mag1)
%title('Gravitational Acceleration Opposed by Solar Powered Thrust rho = 280');
%xlabel('x deviation from L3');ylabel('y deviation')
%figure(6)
%surfc(position_x-L3,position_y,accel_mag2)
%title('Gravitational Acceleration Opposed by Solar Powered Thrust rho = 2800');
%xlabel('x deviation from L3');ylabel('y deviation')
figure(7)
surfc(position_x-L3,position_y,accel_mag3)
title('Gravitational Acceleration Opposed by inefficient Solar Powered Thrust rho = 2800');
xlabel('x deviation from L3');ylabel('y deviation')

figure(8)
surfc(position_x-L3,position_y,accel_test3)
title('Gravitational Acceleration in vicinity of L3');
xlabel('x deviation from L3');ylabel('y deviation')

