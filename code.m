clear all; close all; clc;

%% DH Parameters
syms theta1 theta2 theta3 theta4 theta5 theta6;
%DH = [theta;alpha;r;d] DH is written in this type of matrix with every column representing 4 parameters
DH = [theta1 90 0 330; theta2 0 260 0; theta3+90 90 0 0;...
    theta4 -90 0 290; theta5 90 0 0; theta6 0 0 70];

%% Homogenous matrices
H10 = [cos(theta1) 0 sin(theta1) 0; sin(theta1) 0 -cos(theta1) 0; 0 1 0 330; 0 0 0 1];
H21 = [cos(theta2) -sin(theta2) 0 260*cos(theta2); sin(theta2) cos(theta2) 0 260*sin(theta2); 0 0 1 0; 0 0 0 1];
H32 = [cos(theta3) 0 sin(theta3) 0; sin(theta3) 0 -cos(theta3) 0; 0 1 0 0; 0 0 0 1];
H43 = [cos(theta4) 0 sin(theta4) 0; sin(theta4) 0 -cos(theta4) 0; 0 -1 0 0; 0 0 0 1];
H54 = [cos(theta5) 0 sin(theta5) 0; sin(theta5) 0 -cos(theta5) 0; 0 1 0 0; 0 0 0 1];
H65 = [cos(theta6) -sin(theta6) 0 0; sin(theta6) cos(theta6) 0 0; 0 0 1 70; 0 0 0 1];
H20 = H10*H21;
H30 = H20*H32;
H40 = H30*H43;
H50 = H40*H54;
H60 = H50*H65;

%% Jacobian
A = ([1 0 0;0 1 0;0 0 1]*[0;0;1]);
B = (H60(1:3,4)-[0;0;0]);
v1 = cross(A,B);
C = ([cos(theta1) 0 sin(theta1); sin(theta1) 0 -cos(theta1); 0 1 0])*([0; 0; 1]);
D = (H60(1:3,4)-[0; 0; 330]);
v2 = cross(C,D);
v3 = [cos(theta1)*cos(theta2) -cos(theta1)*sin(theta2) sin(theta1); ...
    cos(theta2)*sin(theta1) -sin(theta1)*sin(theta2) -cos(theta1); ...
    sin(theta2) cos(theta2) 0]*[0;0;1];
E = [cos(theta1)*cos(theta2)*cos(theta3)-cos(theta1)*sin(theta2)*sin(theta3) sin(theta1) cos(theta1)*cos(theta2)*sin(theta3) + cos(theta1)*cos(theta3)*sin(theta2); ...
    cos(theta2)*cos(theta3)*sin(theta1)-sin(theta1)*sin(theta2)*sin(theta3) -cos(theta1) cos(theta2)*sin(theta1)*sin(theta3)+cos(theta3)*sin(theta1)*sin(theta2); ...
    cos(theta2)*sin(theta3)+cos(theta3)*sin(theta2) 0 ...
    sin(theta2)*sin(theta3)-cos(theta2)*cos(theta3)]*[0;0;1];
F = (H60(1:3,4)-[260*cos(theta1)*cos(theta2); 260*cos(theta2)*sin(theta1); 260*sin(theta2) + 330]);
v4 = cross(E,F);
v5 = [H40(1:3,1:3)]*[0;0;1];
G = H50(1:3,1:3)*[0;0;1];
H = H60(1:3,4)-H50(1:3,4);
v6 = cross(G,H);
w1 = [1 0 0; 0 1 0; 0 0 1]*[0;0;1];
w2 = H10(1:3,1:3)*[0;0;1];
w3 = [0;0;0];
w4 = H30(1:3,1:3)*[0;0;1];
w5 = [0;0;0];
w6 = H50(1:3,1:3)*[0;0;1];
jacob = [v1 v2 v3 v4 v5 v6;w1 w2 w3 w4 w5 w6];
jacob1 = simplify(jacob);

%% Inertia matrix
syms m1 m2 m3 m4 m5 m6;
syms x1 x2 x3 x4 x5 x6;
syms y1 y2 y3 y4 y5 y6;
syms z1 z2 z3 z4 z5 z6;

Ixx = (m1*(y1^2+z1^2)) + (m2*(y2^2+z2^2)) + (m3*(y3^2+z3^2)) + ...
    (m4*(y4^2+z4^2)) + (m5*(y5^2+z5^2)) + (m6*(y6^2+z6^2));
Iyy = (m1*(x1^2+z1^2)) + (m2*(x2^2+z2^2)) + (m3*(x3^2+z3^2)) + ...
    (m4*(x4^2+z4^2)) + (m5*(x5^2+z5^2)) + (m6*(x6^2+z6^2));
Izz = (m1*(x1^2+y1^2)) + (m2*(x2^2+y2^2)) + (m3*(x3^2+y3^2)) + ...
    (m4*(x4^2+y4^2)) + (m5*(x5^2+y5^2)) + (m6*(x6^2+y6^2));
Ixy = -(m1*x1*y1) - (m2*x2*y2) - (m1*x3*y3) - (m4*x4*y4) - ...
    (m5*x5*y5) - (m6*x6*y6);
Iyx = -(m1*x1*y1) - (m2*x2*y2) - (m1*x3*y3) - (m4*x4*y4) - ...
    (m5*x5*y5) - (m6*x6*y6);
Ixz = -(m1*x1*z1)-(m2*x2*z2)-(m3*x3*z3)-(m4*x4*z4)...
    -(m5*x5*z5)-(m6*x6*z6);
Izx = -(m1*x1*z1)-(m2*x2*z2)-(m3*x3*z3)-(m4*x4*z4)...
    -(m5*x5*z5)-(m6*x6*z6);
Iyz = -(m1*y1*z1)-(m2*y2*z2)-(m3*y3*z3)-(m4*y4*z4)...
    -(m5*y5*z5)-(m6*y6*z6);
Izy = -(m1*y1*z1)-(m2*y2*z2)-(m3*y3*z3)-(m4*y4*z4)...
    -(m5*y5*z5)-(m6*y6*z6);
IM = [Ixx -Ixy -Ixz; -Ixy Iyy Iyz; -Ixz -Iyz Izz];
Inertia_matirx = simplify(IM);

%% Kinetic Energy
syms Ix1 Ix2 Ix3 Ix4 Ix5 Ix6;
syms L1 L2 L3 L4 L5 L6;

vx1 = -L1*sin(x1)*Ix1;
vx2 = vx1 - L2*sin(x1+x2)*(Ix1+Ix2);
vx3 = vx2 - L3*sin(x1+x2+x3)*(Ix1+Ix2+Ix3);
vx4 = vx3 - L4*sin(x1+x2+x3+x4)*(Ix1+Ix2+Ix3+Ix4);
vx5 = vx4 - L4*sin(x1+x2+x3+x4+x5)*(Ix1+Ix2+Ix3+Ix4+Ix5);
vx6 = vx5 - L5*sin(x1+x2+x3+x4+x5+x6)*(Ix1+Ix2+Ix3+Ix4+Ix5+Ix6);
vy1 = L1*cos(x1)*Ix1;
vy2 = vy1 + L2*cos(x1+x2)*(Ix1+Ix2);
vy3 = vy2 + L3*cos(x1+x2+x3)*(Ix1+Ix2+Ix3);
vy4 = vy3 + L4*cos(x1+x2+x3+x4)*(Ix1+Ix2+Ix3+Ix4);
vy5 = vy4 + L4*cos(x1+x2+x3+x4+x5)*(Ix1+Ix2+Ix3+Ix4+Ix5);
vy6 = vy5 + L5*cos(x1+x2+x3+x4+x5+x6)*(Ix1+Ix2+Ix3+Ix4+Ix5+Ix6);

KE = 0.5*m1*( vx1^2 + vy1^2) + 0.5*m2*( vx2^2 + vy2^2) + 0.5*m3*( vx3^2 + vy4^2) + ...
    0.5*m4*( vx4^2 + vy4^2) + 0.5*m5*( vx5^2 + vy5^2) + 0.5*m6*( vx6^2 + vy6^2);

%% Potential Energy
syms g;
hx1 = L1*cos(x1);
hx2 = hx1+L2*cos(x1+x2);
hx3 = hx2+L3*cos(x1+x2+x3);
hx4 = hx3+L4*cos(x1+x2+x3+x4);
hx5 = hx4+L5*cos(x1+x2+x3+x4+x5);
hx6 = hx5+L6*cos(x1+x2+x3+x4+x5+x6);
hy1 = L1*sin(x1);
hy2 = hy1+L2*sin(x1+x2);
hy3 = hy2+L3*sin(x1+x2+x3);
hy4 = hy3+L4*sin(x1+x2+x3+x4);
hy5 = hy4+L5*sin(x1+x2+x3+x4+x5);
hy6 = hy5+L6*sin(x1+x2+x3+x4+x5+x6);

PE = m1*g*hy1 + m2*g*hy2 + m3*g*hy3 +m4*g*hy4 + m5*g*hy5 + m6*g*hy6;

%% Langrange
Langrange = KE - PE;

%% Differential equation of motion
%sysm q1 q2 q3;
%Torque = Inertia_matrix * [diff(q1,2);diff(q2,2);diff(q3,2)] + ...
%    PE - centrigugal_forces   % i could not figure iiut this part
