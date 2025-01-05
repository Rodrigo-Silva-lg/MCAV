%% State Machine

% Parameters
g_mars = 3.73;
Cd = 1.05;
ro_mars = 0.020;
z_I = [0;0;1];
l_L = 3.5;
w_L = 3.2;
h_L = 2.2;
ml = 1500; 
l_R = 3.0;
w_R = 2.7;
h_R = 2.2;
mr = 1025;
m = 2525;

J_L = [(1/12)*ml*(l_L^2+h_L^2) 0 0  
      0 (1/12)*ml*(w_L^2+h_L^2) 0   
      0 0 (1/12)*ml*(l_L^2+w_L^2)];

J_R = [(1/12)*mr*(l_R^2+h_R^2) 0 0  
      0 (1/12)*mr*(w_R^2+h_R^2) 0
      0 0 (1/12)*mr*(l_R^2+w_R^2)];

cmz = (h_R*mr+h_L*ml)/m;
z = cmz-2.2;

J = [J_L(1,1)+ml*(1.5-z)^2+J_R(1,1)+mr*(cmz-1.1)^2 0 0
     0 J_L(2,2)+ml*(1.5-z)^2+J_R(2,2)+mr*(cmz-1.1)^2 0
     0 0 J_L(3,3)+J_R(3,3)];

beta = -0.5*Cd*ro_mars*l_L*l_L;

%f's & n's
f1 = [0;sin(deg2rad(180-25));cos(deg2rad(180-25))];
p1 = [-1.6;-1.75;0.47525];
n1 = cross(p1,f1);
f2 = [0;sin(deg2rad(180+25));cos(deg2rad(180+25))];
p2 = [-1.6;1.75;0.47525];
n2 = cross(p2,f2);
f3 = [0;sin(deg2rad(180+25));cos(deg2rad(180+25))];
p3 = [1.6;1.75;0.47525];
n3 = cross(p3,f3);
f4 = [0;sin(deg2rad(180-25));cos(deg2rad(180-25))];
p4 = [1.6;-1.75;0.47525];
n4 = cross(p4,f4);

%Model 1
px1 = 0;
py1 = 0;
pz1 = -2100;
p1 = [px1;py1;pz1];
vx1 = 0;
vy1 = 0;
vz1 = 89;
v1 = [vx1;vy1;vz1];
wx1 = 0;  
wy1 = 0;
wz1 = 0;
omega1 = [wx1;wy1;wz1]; 

phi1 = deg2rad(0);
theta1 = deg2rad(0);
psi1 = deg2rad(0);
lbd1 = [phi1;theta1;psi1];
T_1 = 0.2*(m*g_mars+beta*vz1^2);

a_1 = [0, -vx1*sin(theta1)+vz1*cos(theta1), -vx1*sin(psi1)-vy1*cos(psi1)
     -vy1*sin(phi1)-vz1*cos(phi1), 0, vx1*cos(psi1)-vy1*sin(psi1)
      vy1*cos(phi1)-vz1*sin(phi1), vx1*cos(theta1)-vz1*sin(theta1), 0];

b_1 = [(wy1*cos(phi1)-wz1*sin(phi1))*tan(theta1), (wy1*sin(phi1)+wz1*cos(phi1))/((cos(theta1))^2), 0
     -wy1*sin(phi1)-wz1*cos(phi1), 0, 0
      (wy1*cos(phi1)-wz1*sin(phi1))/cos(theta1), ((wy1*sin(phi1)+wz1*cos(phi1))*tan(theta1))/cos(theta1), 0];

A_1 = [ zeros(3), Euler2R(lbd1), a_1, zeros(3)
      zeros(3), [zeros(2,3);0 0 2*beta*vz1/m], skew(g_mars*z_I), skew(v1)
      zeros(3), zeros(3), b_1, Euler2Q(lbd1)
      zeros(3), zeros(3), zeros(3), zeros(3)];

B_1 = [zeros(3,4)
       1/m*f1, 1/m*f2, 1/m*f3, 1/m*f4
       zeros(3,4)
       J^-1*n1, J^-1*n2, J^-1*n3, J^-1*n4];

C = eye(12);

D = zeros(12,4);

%Model 2
px2 = 0;
py2 = 0;
pz2 = -1000;
p2= [px2;py2;pz2];
vx2 = 0;
vy2 = 0;
vz2 = 44.5;
v2 = [vx2;vy2;vz2];
wx2 = 0;  
wy2 = 0;
wz2 = 0;
omega2 = [wx2;wy2;wz2]; 

phi2 = deg2rad(0);
theta2 = deg2rad(0);
psi2 = deg2rad(0);
lbd2 = [phi2;theta2;psi2];
x02 = [p2; v2; lbd2; omega2];
T2 = 0.2*(m*g_mars+beta*vz2^2);

a_2 = [0, -vx2*sin(theta2)+vz2*cos(theta2), -vx2*sin(psi2)-vy2*cos(psi2)
     -vy2*sin(phi2)-vz2*cos(phi2), 0, vx2*cos(psi2)-vy2*sin(psi2)
      vy2*cos(phi2)-vz2*sin(phi2), vx2*cos(theta2)-vz2*sin(theta2), 0];

b_2 = [(wy2*cos(phi2)-wz2*sin(phi2))*tan(theta2), (wy2*sin(phi2)+wz2*cos(phi2))/((cos(theta2))^2), 0
     -wy2*sin(phi2)-wz2*cos(phi2), 0, 0
      (wy2*cos(phi2)-wz2*sin(phi2))/cos(theta2), ((wy2*sin(phi2)+wz2*cos(phi2))*tan(theta2))/cos(theta2), 0];

A_2 = [ zeros(3), Euler2R(lbd2), a_2, zeros(3)
     zeros(3), [zeros(2,3);0 0 2*beta*vz2/m], skew(g_mars*z_I), skew(v2)
     zeros(3), zeros(3), b_2, Euler2Q(lbd2)
     zeros(3), zeros(3), zeros(3), zeros(3)];

B_2 = [zeros(3,4)
   1/m*f1, 1/m*f2, 1/m*f3, 1/m*f4
   zeros(3,4)
   J^-1*n1, J^-1*n2, J^-1*n3, J^-1*n4];

%Model 3
px3 = 0;
py3 = 0;
pz3 = -21;
p3 = [px3;py3;pz3];
vx3 = 0;
vy3 = 0;
vz3 = 0;
v3 = [vx3;vy3;vz3];
wx3 = 0;  
wy3 = 0;
wz3 = 0;
omega3 = [wx3;wy3;wz3]; 

phi3 = deg2rad(0);
theta3 = deg2rad(0);
psi3 = deg2rad(0);
lbd3 = [phi3;theta3;psi3];
x03 = [p3; v3; lbd3; omega3];
T3 = (0.2*m*g_mars);

a_3 = [0, -vx3*sin(theta3)+vz3*cos(theta3), -vx3*sin(psi3)-vy3*cos(psi3)
     -vy3*sin(phi3)-vz3*cos(phi3), 0, vx3*cos(psi3)-vy3*sin(psi3)
      vy3*cos(phi3)-vz3*sin(phi3), vx3*cos(theta3)-vz3*sin(theta3), 0];

b_3 = [(wy3*cos(phi3)-wz3*sin(phi3))*tan(theta3), (wy3*sin(phi3)+wz3*cos(phi3))/((cos(theta3))^2), 0
     -wy3*sin(phi3)-wz3*cos(phi3), 0, 0
      (wy3*cos(phi3)-wz3*sin(phi3))/cos(theta3), ((wy3*sin(phi3)+wz3*cos(phi3))*tan(theta3))/cos(theta3), 0];

A_3 = [ zeros(3), Euler2R(lbd3), a_3, zeros(3)
     zeros(3), [zeros(2,3);0 0 2*beta*vz3/m], skew(g_mars*z_I), skew(v3)
     zeros(3), zeros(3), b_3, Euler2Q(lbd3)
     zeros(3), zeros(3), zeros(3), zeros(3)];

B_3 = [zeros(3,4)
   1/m*f1, 1/m*f2, 1/m*f3, 1/m*f4
   zeros(3,4)
   J^-1*n1, J^-1*n2, J^-1*n3, J^-1*n4];

%Model 4

J = [(1/12)*ml*(l_L^2+h_L^2) 0 0  
      0 (1/12)*ml*(w_L^2+h_L^2) 0   
      0 0 (1/12)*ml*(l_L^2+w_L^2)];

%f's & n's
f1 = [0;sin(deg2rad(180-25));cos(deg2rad(180-25))];
p1 = [-1.6;-1.75;0.47525];
n1 = cross(p1,f1);
f2 = [0;sin(deg2rad(180+25));cos(deg2rad(180+25))];
p2 = [-1.6;1.75;0.47525];
n2 = cross(p2,f2);
f3 = [0;sin(deg2rad(180+25));cos(deg2rad(180+25))];
p3 = [1.6;1.75;0.47525];
n3 = cross(p3,f3);
f4 = [0;sin(deg2rad(180-25));cos(deg2rad(180-25))];
p4 = [1.6;-1.75;0.47525];
n4 = cross(p4,f4);

T4 = (ml*g_mars*0.2);
px4 = 0;
py4 = 0;
pz4 = -21;
p4 = [px4;py4;pz4];
vx4 = 0;
vy4 = 10;
vz4 = 0;
v4 = [vx4;vy4;vz4];
wx4 = 0;  
wy4 = 0;
wz4 = 0;
omega4 = [wx4;wy4;wz4]; 
phi4 = deg2rad(1);
theta4 = deg2rad(0);
psi4 = deg2rad(0);
lbd4 = [phi4;theta4;psi4];
x04 = [p4; v4; lbd4; omega4];

a_4 = [0, -vx4*sin(theta4)+vz4*cos(theta4), -vx4*sin(psi4)-vy4*cos(psi4)
     -vy4*sin(phi4)-vz4*cos(phi4), 0, vx4*cos(psi4)-vy4*sin(psi4)
      vy4*cos(phi4)-vz4*sin(phi4), vx4*cos(theta4)-vz4*sin(theta4), 0];

b_4 = [(wy4*cos(phi4)-wz4*sin(phi4))*tan(theta4), (wy4*sin(phi4)+wz4*cos(phi4))/((cos(theta4))^2), 0
     -wy4*sin(phi4)-wz4*cos(phi4), 0, 0
      (wy4*cos(phi4)-wz4*sin(phi4))/cos(theta4), ((wy4*sin(phi4)+wz4*cos(phi4))*tan(theta4))/cos(theta4), 0];

c_4 = [0, -g_mars, 0
     g_mars*cos(phi4), 0, 0
     -g_mars*sin(phi4), 0, 0];


A_4 = [ zeros(3), Euler2R(lbd4), a_4, zeros(3)
     zeros(3), zeros(3), c_4, skew(v4)
     zeros(3), zeros(3), b_4, Euler2Q(lbd4)
     zeros(3), zeros(3), zeros(3), zeros(3)];

B_4 = [zeros(3,4)
   1/m*f1, 1/m*f2, 1/m*f3, 1/m*f4
   zeros(3,4)
   J^-1*n1, J^-1*n2, J^-1*n3, J^-1*n4];

%Controllers used
%H_inf for situation 1
B1_1 = zeros(12);
B2_1 = B_1;
W1_1 = sqrt(Q1);
W2_1 = sqrt(R);
C1_1 = [W1_1];
D11_1 = zeros(12);
D12_1 = [zeros(8,4);W2_1];
C2_1 = -eye(12);
D21_1 = eye(12);
D22_1 = zeros(12,4);
A01 = A_1;
B01 = [B1_1,B2_1];
C01 = [C1_1;C2_1];
D01 = [D11_1 D12_1; D21_1 D22_1];

P1 = ss(A01,B01,C01,D01);

%Hinf for situation 2
B1_2 = zeros(12);
B2_2 = B_2;
W1_2 = sqrt(Q2);
W2_2 = sqrt(R);
C1_2 = [W1_2];
D11_2 = zeros(12);
D12_2 = [zeros(8,4);W2_2];
C2_2 = -eye(12);
D21_2 = eye(12);
D22_2 = zeros(12,4);
A02 = A_2;
B02 = [B1_2,B2_2];
C02 = [C1_2;C2_2];
D02 = [D11_2 D12_2; D21_2 D22_2];

P2 = ss(A02,B02,C02,D02);

%LQR for situation 3 (hover)
Q = blkdiag(1,1,10,1,1,100,0.01*eye(3),0.001*eye(3));
R = 0.1*eye(4);
Klqr3 = lqr(A_3,B_3,Q,R);
lbd_CL_lqr3 = eig(A_3-B_3*Klqr3);

%Hinf for situation 4
B1_4 = zeros(12);
B2_4 = B_4;
W1_4 = sqrt(Q4);
W2_4 = sqrt(R);
C1_4 = [W1_4];
D11_4 = zeros(12);
D12_4 = [zeros(8,4);W2_4];
C2_4 = -eye(12);
D21_4 = eye(12);
D22_4 = zeros(12,4);
A04 = A_4;
B04 = [B1_4,B2_4];
C04 = [C1_4;C2_4];
D04 = [D11_4 D12_4; D21_4 D22_4];

P4 = ss(A04,B04,C04,D04);

nmeas = 12;
ncont = 4;
% tests for each case/controller:
[Kinf1,CLinf1,gammainf1,info_inf1] = hinfsyn(P1,nmeas,ncont);
poles_CLinf1 = pole(CLinf1);
if any(real(poles_CLinf1) >= 0), disp('CL system with Hinf1 controller not stable'); else, disp('CL system with Hinf1 controller is stable'); end

[Kinf2,CLinf2,gammainf2,info_inf2] = hinfsyn(P2,nmeas,ncont);
poles_CLinf2 = pole(CLinf2);
if any(real(poles_CLinf2) >= 0), disp('CL system with Hinf2 controller not stable'); else, disp('CL system with Hinf2 controller is stable'); end

if any(real(lbd_CL_lqr3) >= 0), disp('CL system with LQR3 not stable'); else, disp('CL system with LQR3 is stable'); end

[Kinf4,CLinf4,gammainf4,info_inf4] = hinfsyn(P4,nmeas,ncont);
poles_CLinf4 = pole(CLinf4);
if any(real(poles_CLinf4) >= 0), disp('CL system with Hinf4 controller not stable'); else, disp('CL system with Hinf4 controller is stable'); end

% State Machine - Implementation
Dt = 0.01;
t = 0:Dt:120;
NSim = length(t);
nx = 12;
nu = 4;
r = zeros(nx,NSim);
x = zeros(nx,NSim);
u = zeros(nu,NSim);
i_ctr = 0;
aux = 10000000;
x(:,1) = [0;0;-2100;0;0;89;0;0;0;0;0;0];

for k = 1:NSim
    if x(3,k) >= -2100 && x(3,k) < -1000 && k<=aux
        i_ctr = 1;
        x(:,1) = x1;
        r(:,k) = [0;0;-1000;0;0;0;0;0;0;0;0;0]*(t(k)>=0);    

    elseif x(3,k) >= -1000 && x(3,k) < -21 && k<=aux
        i_ctr = 2;
        x(:,1)=x2;
        r(:,k) = [0;0;-21;0;0;0;0;0;0;0;0;0]*(t(k)>=0);      

    elseif x(3,k) >= -21 && k<=aux
        i_ctr = 3;
        x(:,1)=x3;
        r(:,k) = [0;0;-21;0;0;0;0;0;0;0;0;0]*(t(k)>=0);    
        if aux==10000000
            aux = k + 1499;
        end

    elseif k>=aux
        i_ctr = 4;
        x(:,1)=x4;
        r(:,k) = [0;500;-21;0;0;0;0;0;0;0;0;0]*(t(k)>=0);     

    end

    y(:, k) = C * x(:, k);

    switch i_ctr
        case 1 
            u(:, k) = -Kinf1.C * [r(:, k) - y(:, k)];

            x_dot = A_1*x(:,k)+B_1*u(:,k);
            xp = x(:,k)+Dt*x_dot;
            if k < NSim
                x(:,k+1) = xp;
            end
        case 2
            u(:, k) = -Kinf2.C * (r(:, k) - y(:, k));

            x_dot = A_2 * x(:, k) + B_2 * u(:, k); 
            xp = x(:,k) + Dt * x_dot;
            if k < NSim
                x(:,k + 1) = xp;
            end
        case 3 
            u(:,k) = -Klqr3*(y(:,k) - r(:,k));

            x_dot = A_3 * x(:,k) + B_3 * u(:,k);
            xp = x(:,k) + Dt * x_dot; 
            if k < NSim
                x(:, k + 1) = xp;
            end
        case 4
            u(:, k) = -Kinf4.C * (r(:,k) - y(:,k));

            x_dot = A_4 * x(:,k) + B_4 * u(:,k); 
            xp = x(:,k) + Dt * x_dot; 
            if k < NSim
                x(:, k + 1) = xp;
            end
    end

 end

figure(1);
plot(t,y(2,:),t,y(3,:));
ylabel('p(t)');
xlabel('Tempo [s]');
legend('p_y','p_z');
grid on;

figure(2);
plot(t,u(1,:),t,u(2,:),t,u(3,:),t,u(4,:));
ylabel('u(t)');
xlabel('t [s]');
legend('u_1', 'u_2', 'u_3', 'u_4');
grid on;        

%% Auxiliary Functions

function rot = Euler2R(ang)


rotx = [	1	,	0			,	 0
			0	,	cos(ang(1))	,	-sin(ang(1))
			0	,	sin(ang(1))	,	 cos(ang(1))	];
roty = [	 cos(ang(2))	,	0	,	sin(ang(2))
			0				,	1	,	0
			-sin(ang(2))	,	0	,	cos(ang(2))	];
rotz = [	cos(ang(3))	,	-sin(ang(3))	,	0
			sin(ang(3))	,	 cos(ang(3))	,	0
			0			,	 0				,	1	];
rot = rotz*roty*rotx;

end

function Q = Euler2Q(l)

phi = l(1); theta = l(2);
Q = [  1 sin(phi)*tan(theta) cos(phi)*tan(theta);
       0 cos(phi)           -sin(phi);
       0 sin(phi)/cos(theta) cos(phi)/cos(theta)];
end

function X=skew(x)


n = length(x);
if n == 3
    X = [   0      -x(3)    x(2)
            x(3)    0      -x(1)
           -x(2)    x(1)    0     ];
elseif n == 1
    X = [   0      -x(1)
            x(1)    0    ];
else
    error('SKEW function not implemented for input dimensions other than 1 or 3 (i.e., so(2) and so(3)).');
end

end 