%% 1 LANDER NONLINEAR MODEL


%Parameters

g_mars = 3.73;
Cd = 1.05;
ro_mars=0.020;

xr = 3;
yr = 2.7;
zr = 2.2;
mr = 1025;

jr = [(1/12)*mr*(yr^2+zr^2) 0 0
      0 (1/12)*mr*(xr^2+zr^2) 0
      0 0 (1/12)*mr*(yr^2+xr^2)];

xl=3.5;
yl=3.2;
zl = 3;
ml = 1500;

jl = [(1/12)*ml*(yl^2+zl^2) 0 0
      0 (1/12)*ml*(xl^2+zl^2) 0
      0 0 (1/12)*ml*(yl^2+xl^2)];

m = mr+ml;

cmz = (zr*mr+zl*ml)/m;

z = cmz-2.2;

J = [jl(1,1)+ml*(1.5-z)^2+jr(1,1)+mr*(cmz-1.1)^2 0 0
     0 jl(2,2)+ml*(1.5-z)^2+jr(2,2)+mr*(cmz-1.1)^2 0
     0 0 jl(3,3)+jr(3,3)];

%f's & n's assuming T=1

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

%Simulation

zI=[0;0;1];
Dt=0.01;
t=0:Dt:60;
Nsim=length(t);
x0=zeros(12,1);
x=zeros(12,Nsim);
y=zeros(4,Nsim);
x(:,1) = x0;

T1 = 0.2*m*g_mars;
T2 = 0.2*m*g_mars;
T3 = 0.2*m*g_mars;
T4 = 0.2*m*g_mars;

for k = 1:Nsim

  p=x(1:3,k);
  v=x(4:6,k);
  lambda=x(7:9,k);
  omega=x(10:12,k);
  R=Euler2R(lambda);
  Q=Euler2Q(lambda);
  T=T1*f1+T2*f2+T3*f3+T4*f4;
  np=T1*n1+T2*n2+T3*n3+T4*n4;

  p_dot=R*v;
  lambda_dot=Q*omega;
  v_dot=-skew(omega)*v+(g_mars*R.'*zI)-(0.5*Cd*ro_mars*xr*yr*zI*v(3)^2/m)+(T/m);
  omega_dot=(-J^-1*skew(omega)*J*omega)+(J^-1*np);
  x_dot=[p_dot;v_dot;lambda_dot;omega_dot];


  x(:,k+1)=x(:,k)+Dt*x_dot;
  y(1,k)=x(1,k); %px
  y(2,k)=x(2,k); %py
  y(3,k)=x(3,k); %pz
  y(4,k)=x(6,k); %vz

end

figure(33);
subplot(3,1,1)
plot(t,y(1,:));
legend('px')
axis([0 60 -1 1])
grid on
subplot(3,1,2)
plot(t,y(2,:));
legend('py');
axis([0 60 -1 1])
grid on
subplot(3,1,3)
plot(t,y(3,:));
legend('pz');
grid on

figure(1);
plot(1)
plot(t,y(4,:));
legend('vz');
grid on


% 2 LANDER LINEARIZED MODEL

%% MODEL 1

%Parameters
Dt = 0.01;
t = 0:Dt:60;
Nsim = length(t);
T = (m*g_mars*0.25);
u = ones(Nsim,1)*[T,T,T,T];
beta = -0.5*Cd*ro_mars*yl*xl;

%Equilibrium Conditions
px = 0;
py = 0;
pz = 0;
p = [px;py;pz];
wx = 0;  
wy = 0;
wz = 0;
omega = [wx;wy;wz]; 
vx = 0;
vy = 0;
vz = 89;
v = [vx;vy;vz];
phi = deg2rad(0);
theta = deg2rad(0);
psi = deg2rad(0);
lbd = [phi;theta;psi];
x0 = [p; v; lbd; omega];
T = 0.25*(m*g_mars+beta*vz^2);
u = ones(Nsim,1)*[T,T,T,T];


a = [0, -vx*sin(theta)+vz*cos(theta), -vx*sin(psi)-vy*cos(psi)
     -vy*sin(phi)-vz*cos(phi), 0, vx*cos(psi)-vy*sin(psi)
      vy*cos(phi)-vz*sin(phi), vx*cos(theta)-vz*sin(theta), 0];

b = [(wy*cos(phi)-wz*sin(phi))*tan(theta), (wy*sin(phi)+wz*cos(phi))/((cos(theta))^2), 0
     -wy*sin(phi)-wz*cos(phi), 0, 0
      (wy*cos(phi)-wz*sin(phi))/cos(theta), ((wy*sin(phi)+wz*cos(phi))*tan(theta))/cos(theta), 0];

A = [ zeros(3), Euler2R(lbd), a, zeros(3)
     zeros(3), [zeros(2,3);0 0 2*beta*vz/m], skew(g_mars*zI), skew(v)
     zeros(3), zeros(3), b, Euler2Q(lbd)
     zeros(3), zeros(3), zeros(3), zeros(3)];

B=[zeros(3,4)
   1/m*f1, 1/m*f2, 1/m*f3, 1/m*f4
   zeros(3,4)
   J^-1*n1, J^-1*n2, J^-1*n3, J^-1*n4];

C = eye(12);

D = zeros(12,4);

sys = ss(A,B,C,D);
y_L = lsim(sys,u,t,x0);

% Controllability, observability and stability
[M,J] = jordan(A),
if rank(ctrb(A,B)) < size(A,1), disp('The system is not controllable');
else 
    disp('The system is controllable');
end
if rank(obsv(A,C)) < size(A,1), disp('The system is not observable');
else 
    disp('The system is observable');
end


figure(6);
plot(t,y_L(:,1),'k',t,y_L(:,2),'r',t,y_L(:,3),'g',t,y_L(:,4),'b');
legend('p_x','p_y','p_z','v_z')
grid on;

% SVD Frequency Analysis
figure(2);
sigma(sys);
title('Sigma Values (Vz = 89 m/s, Vx = Vy = 0)');
grid on;

% Checking stability with LQR controller
Q = blkdiag(1,10,100,1,100,10,0.01*eye(3),0.001*eye(3));
R = 0.1*eye(4);
Klqr = lqr(A,B,Q,R);
lbd_CL_lqr = eig(A-B*Klqr); 
X = A-B*Klqr;
CLlqr = ss(X,B,C,D);
if any(real(lbd_CL_lqr) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end

% LQR Controller - Implementation
Dt = 0.01;
t = 0:Dt:60;
NSim = length(t);
nx = 12;
nu = 4;
x = zeros(nx, NSim);
u = zeros(nu, NSim);
x(:, 1) = [0; 0; -2100; 0; 0; 89; 0; 0; 0; 0; 0; 0]; % Initial Conditions

% Simulation loop with LQR controller
for k = 1:NSim
    % Measurements
    y1(:, k) = C * x(:, k);
    
    % Control using LQR
    u(:, k) = -Klqr * x(:, k);
    
    % System Dynamics
    x_dot = A * x(:, k) + B * u(:, k);
    xp = x(:, k) + Dt * x_dot;
    if k < NSim
        x(:, k + 1) = xp;
    end
end

figure(30);
plot(t, u);
grid on;
xlabel('t [s]');
ylabel('u(t)');
legend('u_1', 'u_2', 'u_3', 'u_4');
title('Control Inputs with LQR Controller');

figure(31);
plot(t, y1(3, :)); % p_z
grid on;
xlabel('t [s]');
ylabel('p_z');
legend('p_z');
title('Position p_z with LQR Controller');


% H_inf Controller - Definition

A0 = A;
B1 = zeros(12,12);
B2 = B;
W1 = sqrt(Q);
W2 = sqrt(R);
C1 = [W1];
D11 = zeros(12,12);
D12 = [zeros(8,4); W2];
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);

B0 = [B1, B2];
C0 = [C1; C2];
D0 = [D11, D12; D21, D22];

P = ss(A0, B0, C0, D0);

nmeas = 12; 
ncont = 4; 

% tests on P:
if (12 - rank(ctrb(A0,B2))) > 0, disp('A1.1 on P: system uncontrolable'); else disp('A1.1 on P: OK'); end
if (12 - rank(obsv(A0,C2))) > 0, disp('A1.2 on P: system unobservable'); else disp('A1.2 on P: OK'); end
if (size(D12,2) - rank(D12)) > 0, disp('A2.1 on P: D12 is column rank deficient'); else disp('A2.1 on P: OK'); end
if (size(D21,1) - rank(D21)) > 0, disp('A2.2 on P: D21 is row rank deficient'); else disp('A2.1 on P: OK'); end
syms w real; 
Aux1 = [A0 - j*w*eye(size(A0)) , B2 ; C1 , D12];
if (size(Aux1,2) - rank(Aux1)) > 0,  disp('A3 on P: matrix is column rank deficient'); else disp('A3 on P: OK'); end
Aux2 = [A0 - j*w*eye(size(A0)) , B1 ; C2 , D21];
if (size(Aux2,1) - rank(Aux2)) > 0,  disp('A4 on P: matrix is column rank deficient'); else disp('A4 on P: OK'); end

[Kinf,CLinf,gammainf,info_inf] = hinfsyn(P,nmeas,ncont);
poles_CLinf = pole(CLinf);
if any(real(poles_CLinf) >= 0), disp('CL system with Hinf controller not stable'); else, disp('CL system with Hinf controller is stable'); end

% Hinf controller - Implementation
Dt = 0.01;
t = 0:Dt:60;
NSim = length(t);
nx = 12;
nu = 4;
x = zeros(nx,NSim);
xu = zeros(nx,NSim);
u = zeros(nu,NSim);
x(:,1) = [0;0;-2100;0;0;89;0;0;0;0;0;0]; %Initial Conditions
C=eye(12);
for k = 1:NSim

    % Measurements
    y1(:,k) = C*x(:,k);
    % Control
    v = -[-y1(:,k)];
    u(:,k) = Kinf.C*v; % approximation for LQR-like performance of Hinf
    
    % Simulation:
    x_dot = A*x(:,k) + B*u(:,k); 
    xp = x(:,k) + Dt*x_dot; 
    if k < NSim
        x(:,k+1) = xp;
    end
end

figure(20);
plot(t, u);
grid on;
xlabel('t [s]');
ylabel('u(t)');
legend('u_1', 'u_2', 'u_3', 'u_4');
title('Control Inputs with Hinf Controller');

figure(21);
plot(t, y1(3, :));
grid on;
xlabel('t [s]');
ylabel('p_z');
legend('p_z');
title('Position p_z with Hinf Controller');

figure(110);
sigmaplot(CLlqr,'k',Kinf);
grid on;
title('Sigma Values - CL (modelo 1)');
legend('LQR','H\infty');

%% MODEL 2

%Parameters
Dt = 0.01;
t = 0:Dt:60;
Nsim = length(t);
T = (m*g_mars*0.25);
u = ones(Nsim,1)*[T,T,T,T];
beta = -0.5*Cd*ro_mars*yl*xl;

J = [jl(1,1)+ml*(1.5-z)^2+jr(1,1)+mr*(cmz-1.1)^2 0 0
     0 jl(2,2)+ml*(1.5-z)^2+jr(2,2)+mr*(cmz-1.1)^2 0
     0 0 jl(3,3)+jr(3,3)];

%Equilibrium Conditions
px = 0;
py = 0;
pz = 0;
p = [px;py;pz];
wx = 0;  
wy = 0;
wz = 0;
omega = [wx;wy;wz]; 
vx = 0;
vy = 0;
vz = 44.5;
v = [vx;vy;vz];
phi = deg2rad(0);
theta = deg2rad(0);
psi = deg2rad(0);
lbd = [phi;theta;psi];
x0 = [p; v; lbd; omega];
T = 0.25*(m*g_mars+beta*vz^2);
u = ones(Nsim,1)*[T,T,T,T];


a = [0, -vx*sin(theta)+vz*cos(theta), -vx*sin(psi)-vy*cos(psi)
     -vy*sin(phi)-vz*cos(phi), 0, vx*cos(psi)-vy*sin(psi)
      vy*cos(phi)-vz*sin(phi), vx*cos(theta)-vz*sin(theta), 0];

b = [(wy*cos(phi)-wz*sin(phi))*tan(theta), (wy*sin(phi)+wz*cos(phi))/((cos(theta))^2), 0
     -wy*sin(phi)-wz*cos(phi), 0, 0
      (wy*cos(phi)-wz*sin(phi))/cos(theta), ((wy*sin(phi)+wz*cos(phi))*tan(theta))/cos(theta), 0];

A = [ zeros(3), Euler2R(lbd), a, zeros(3)
     zeros(3), [zeros(2,3);0 0 2*beta*vz/m], skew(g_mars*zI), skew(v)
     zeros(3), zeros(3), b, Euler2Q(lbd)
     zeros(3), zeros(3), zeros(3), zeros(3)];

B=[zeros(3,4)
   1/m*f1, 1/m*f2, 1/m*f3, 1/m*f4
   zeros(3,4)
   J^-1*n1, J^-1*n2, J^-1*n3, J^-1*n4];

C = eye(12);

D = zeros(12,4);

sys = ss(A,B,C,D);
y_L = lsim(sys,u,t,x0);


% Controllability, observability and stability
[M,J] = jordan(A),
if rank(ctrb(A,B)) < size(A,1), disp('The system is not controllable');
else 
    disp('The system is controllable');
end
if rank(obsv(A,C)) < size(A,1), disp('The system is not observable');
else 
    disp('The system is observable');
end


figure(7);
plot(t,y_L(:,1),'k',t,y_L(:,2),'r',t,y_L(:,3),'g',t,y_L(:,4),'b');
legend('p_x','p_y','p_z','v_z')
grid on;
xlabel('t [s]');
% SVD Frequency Analysis
figure(3);
sigma(sys);
title('Sigma Values (Vz = 44.5 m/s, Vx = Vy = 0)');
grid on;

% Checking stability with LQR controller
Q = blkdiag(1,10,100,1,100,10,0.01*eye(3),0.001*eye(3));
R = 0.1*eye(4);
Klqr = lqr(A,B,Q,R);
lbd_CL_lqr = eig(A-B*Klqr); 
X = A-B*Klqr;
CLlqr = ss(X,B,C,D);
if any(real(lbd_CL_lqr) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end

% LQR Controller - Implementation
Dt = 0.01;
t = 0:Dt:60;
NSim = length(t);
nx = 12;
nu = 4;
x = zeros(nx, NSim);
u = zeros(nu, NSim);
x(:, 1) = [0; 0; -1000; 0; 0; 44.5; 0; 0; 0; 0; 0; 0]; % Initial Conditions
C=eye(12);

% Simulation loop with LQR controller
for k = 1:NSim
    % Measurements
    y1(:, k) = C * x(:, k);
    
    % Control using LQR
    u(:, k) = -Klqr * x(:, k);
    
    % System Dynamics
    x_dot = A * x(:, k) + B * u(:, k);
    xp = x(:, k) + Dt * x_dot;
    if k < NSim
        x(:, k + 1) = xp;
    end
end

figure(37);
plot(t, u);
grid on;
xlabel('t [s]');
ylabel('u(t)');
legend('u_1', 'u_2', 'u_3', 'u_4');
title('Control Inputs with LQR Controller');

figure(38);
plot(t, y1(3, :));
grid on;
xlabel('t [s]');
ylabel('p_z');
legend('p_z');
title('Position p_z with LQR Controller');

% H_inf Controller - Definition

A0 = A;
B1 = zeros(12,12);
B2 = B;
W1 = sqrt(Q);
W2 = sqrt(R);
C1 = [W1];
D11 = zeros(12,12);
D12 = [zeros(8,4); W2];
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);

B0 = [B1, B2];
C0 = [C1; C2];
D0 = [D11, D12; D21, D22];

P = ss(A0, B0, C0, D0);

nmeas = 12; 
ncont = 4; 

% tests on P:
if (12 - rank(ctrb(A0,B2))) > 0, disp('A1.1 on P: system uncontrolable'); else disp('A1.1 on P: OK'); end
if (12 - rank(obsv(A0,C2))) > 0, disp('A1.2 on P: system unobservable'); else disp('A1.2 on P: OK'); end
if (size(D12,2) - rank(D12)) > 0, disp('A2.1 on P: D12 is column rank deficient'); else disp('A2.1 on P: OK'); end
if (size(D21,1) - rank(D21)) > 0, disp('A2.2 on P: D21 is row rank deficient'); else disp('A2.1 on P: OK'); end
syms w real; 
Aux1 = [A0 - j*w*eye(size(A0)) , B2 ; C1 , D12];
if (size(Aux1,2) - rank(Aux1)) > 0,  disp('A3 on P: matrix is column rank deficient'); else disp('A3 on P: OK'); end
Aux2 = [A0 - j*w*eye(size(A0)) , B1 ; C2 , D21];
if (size(Aux2,1) - rank(Aux2)) > 0,  disp('A4 on P: matrix is column rank deficient'); else disp('A4 on P: OK'); end

[Kinf,CLinf,gammainf,info_inf] = hinfsyn(P,nmeas,ncont);
poles_CLinf = pole(CLinf);
if any(real(poles_CLinf) >= 0), disp('CL system with Hinf controller not stable'); else, disp('CL system with Hinf controller is stable'); end

% Hinf controller - Implementation
Dt = 0.01;
t = 0:Dt:60;
NSim = length(t);
nx = 12;
nu = 4;
x = zeros(nx,NSim);
xu = zeros(nx,NSim);
u = zeros(nu,NSim);
x(:,1) = [0;0;-1000;0;0;44.5;0;0;0;0;0;0]; %Initial Conditions
C=eye(12);
for k = 1:NSim

    % Measurements
    y1(:,k) = C*x(:,k);
    % Control
    v = [y1(:,k)];
    u(:,k) = Kinf.C*v; 
    
    % Simulation
    x_dot = A*x(:,k) + B*u(:,k); 
    xp = x(:,k) + Dt*x_dot; 
    if k < NSim
        x(:,k+1) = xp;
    end
end

figure(22);
plot(t,u);
grid on;
xlabel('t [s]');
ylabel('u(t)');
legend('u_1', 'u_2', 'u_3', 'u_4');
title('Control Inputs with Hinf Controller');

figure(23);
plot(t,y1(3,:));
grid on;
xlabel('t [s]');
ylabel('p_z');
legend('p_z');
title('Position p_z with Hinf Controller');

figure(111);
sigmaplot(CLlqr,'k',Kinf);
grid on;
title('Sigma Values - CL (modelo 2)');
legend('LQR','H\infty');

%% MODEL 3

%Parameters
Dt = 0.01;
t = 0:Dt:60;
Nsim = length(t);
T = (m*g_mars*0.25);
u = ones(Nsim,1)*[T,T,T,T];
beta = -0.5*Cd*ro_mars*yl*xl;

%Equilibrium Conditions
px = 0;
py = 0;
pz = 0;
p = [px;py;pz];
wx = 0;  
wy = 0;
wz = 0;
omega = [wx;wy;wz]; 
vx = 0;
vy = 0;
vz = 0;
v = [vx;vy;vz];
phi = deg2rad(0);
theta = deg2rad(0);
psi = deg2rad(0);
lbd = [phi;theta;psi];
x0 = [p; v; lbd; omega];

a = [0, -vx*sin(theta)+vz*cos(theta), -vx*sin(psi)-vy*cos(psi)
     -vy*sin(phi)-vz*cos(phi), 0, vx*cos(psi)-vy*sin(psi)
      vy*cos(phi)-vz*sin(phi), vx*cos(theta)-vz*sin(theta), 0];

b = [(wy*cos(phi)-wz*sin(phi))*tan(theta), (wy*sin(phi)+wz*cos(phi))/((cos(theta))^2), 0
     -wy*sin(phi)-wz*cos(phi), 0, 0
      (wy*cos(phi)-wz*sin(phi))/cos(theta), ((wy*sin(phi)+wz*cos(phi))*tan(theta))/cos(theta), 0];

c = [0, -g_mars, 0
     g_mars*cos(phi), 0, 0
     -g_mars*sin(phi), 0, 0];


A = [ zeros(3), Euler2R(lbd), a, zeros(3)
     zeros(3), zeros(3), c, skew(v)
     zeros(3), zeros(3), b, Euler2Q(lbd)
     zeros(3), zeros(3), zeros(3), zeros(3)];

B=[zeros(3,4)
   1/m*f1, 1/m*f2, 1/m*f3, 1/m*f4
   zeros(3,4)
   jl^-1*n1, jl^-1*n2, jl^-1*n3, jl^-1*n4];

C = eye(12);

D = zeros(12,4);

sys = ss(A,B,C,D);
y_L = lsim(sys,u,t,x0);

% Controllability, observability and stability
[M,J] = jordan(A),
if rank(ctrb(A,B)) < size(A,1), disp('The system is not controllable');
else 
    disp('The system is controllable');
end
if rank(obsv(A,C)) < size(A,1), disp('The system is not observable');
else 
    disp('The system is observable');
end


figure(8);
plot(t,y_L(:,1),'k',t,y_L(:,2),'r',t,y_L(:,3),'g',t,y_L(:,4),'y',t,y_L(:,5),'b',t,y_L(:,6),'c');
legend('p_x','p_y','p_z','v_x','v_y','v_z')
grid on;
xlabel('t [s]');

% SVD Frequency Analysis
figure(4);
sigma(sys);
title('Sigma Values (Vz = Vx = Vy = 0)');
grid on;

% Checking stability with LQR controller
Q = blkdiag(1,10,100,1,100,10,0.01*eye(3),0.001*eye(3));
R = 0.1*eye(4);
Klqr = lqr(A,B,Q,R);
lbd_CL_lqr = eig(A-B*Klqr); 
X = A-B*Klqr;
CLlqr = ss(X,B,C,D);
if any(real(lbd_CL_lqr) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end

% LQR Controller - Implementation
Dt = 0.01;
t = 0:Dt:60;
NSim = length(t);
x_LQR = zeros(nx, NSim); 
u_LQR = zeros(nu, NSim); 
x_LQR(:,1) = [0; 0; -21; 0; 0; 0; 0; 0; 0; 0; 0; 0]; 
r_LQR = [0; 0; -21; 0; 0; 0; 0; 0; 0; 0; 0; 0] * ones(1, NSim);

% Simulation - LQR controller
for k = 1:NSim
    % Control:
    u_LQR(:,k) = -Klqr * (x_LQR(:,k) - r_LQR(:,k)); 
    
    % Simulation:
    x_dot_LQR = A * x_LQR(:,k) + B * u_LQR(:,k); 
    x_next_LQR = x_LQR(:,k) + Dt * x_dot_LQR; 
    if k < NSim
        x_LQR(:,k+1) = x_next_LQR; 
    end
end

% Results - LQR controller
figure(80);
plot(t, u_LQR);
grid on;
xlabel('t [s]');
ylabel('u(t)');
legend('u_1', 'u_2', 'u_3', 'u_4');
title('Control Inputs with LQR Controller');

figure(81);
plot(t, x_LQR(2,:)); 
grid on;
xlabel('t [s]');
ylabel('\theta (roll)');
legend({'\theta (roll)'});
title('Roll Angle (theta) with LQR Controller');

% H_inf Controller - Definition

A0 = A;
B1 = zeros(12,12);
B2 = B;
W1 = sqrt(Q);
W2 = sqrt(R);
C1 = [W1];
D11 = zeros(12,12);
D12 = [zeros(8,4); W2];
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);

B0 = [B1, B2];
C0 = [C1; C2];
D0 = [D11, D12; D21, D22];

P = ss(A0, B0, C0, D0);

nmeas = 12; 
ncont = 4; 

% tests on P:
if (12 - rank(ctrb(A0,B2))) > 0, disp('A1.1 on P: system uncontrolable'); else disp('A1.1 on P: OK'); end
if (12 - rank(obsv(A0,C2))) > 0, disp('A1.2 on P: system unobservable'); else disp('A1.2 on P: OK'); end
if (size(D12,2) - rank(D12)) > 0, disp('A2.1 on P: D12 is column rank deficient'); else disp('A2.1 on P: OK'); end
if (size(D21,1) - rank(D21)) > 0, disp('A2.2 on P: D21 is row rank deficient'); else disp('A2.1 on P: OK'); end
syms w real; 
Aux1 = [A0 - j*w*eye(size(A0)) , B2 ; C1 , D12];
if (size(Aux1,2) - rank(Aux1)) > 0,  disp('A3 on P: matrix is column rank deficient'); else disp('A3 on P: OK'); end
Aux2 = [A0 - j*w*eye(size(A0)) , B1 ; C2 , D21];
if (size(Aux2,1) - rank(Aux2)) > 0,  disp('A4 on P: matrix is column rank deficient'); else disp('A4 on P: OK'); end

[Kinf,CLinf,gammainf,info_inf] = hinfsyn(P,nmeas,ncont);
poles_CLinf = pole(CLinf);
if any(real(poles_CLinf) >= 0), disp('CL system with Hinf controller not stable'); else, disp('CL system with Hinf controller is stable'); end

% Hinf controller - Implementation
Dt = 0.01;
t = 0:Dt:60;
r = [0;0;-21;0;0;0;0;0;0;0;0;0]*(t>=0);
NSim = length(t);
nx = 12;
nu = 4;
x = zeros(nx,NSim);
xu = zeros(nx,NSim);
u = zeros(nu,NSim);
x(:,1) = [0;0;-21;0;0;0;0;0;0;0;0;0];
C=eye(12);

r = [0;0;-21;0;0;0;0;0;0;0;0;0] * ones(1, length(t));
y1 = zeros(size(C, 1), length(t)); 

for k = 1:NSim
    % Get measurements:
    y1(:,k) = C * x(:,k);

    % Get control action:
    v = -[(r(:,k) - y1(:,k))];
    u(:,k) = Kinf.C * v; 
    
    % Simulate system:
    x_dot = A * x(:,k) + B * u(:,k);
    xp = x(:,k) + Dt * x_dot; 
    if k < NSim
        x(:,k+1) = xp;
    end
end

% Plotting results
figure(24);
plot(t, u);
grid on;
xlabel('t [s]');
ylabel('u(t)');
legend('u_1', 'u_2', 'u_3', 'u_4');
title('Control Inputs with Hinf Controller');

figure(25);
plot(t, y1(2,:)); 
grid on;
xlabel('t [s]');
ylabel('\theta (roll)');
legend({'\theta (roll)'});
title('Roll Angle (theta) with Hinf Controller');

figure(112);
sigmaplot(CLlqr,'k',Kinf);
grid on;
title('Sigma Values - CL (modelo 3)');
legend('LQR','H\infty');

%% MODEL 4

%Parameters
Dt = 0.01;
t = 0:Dt:60;
Nsim = length(t);
T = (m*g_mars*0.25);
u = ones(Nsim,1)*[T,T,T,T];
beta = -0.5*Cd*ro_mars*yl*xl;

%Equilibrium Conditions
px = 0;
py = 0;
pz = 0;
p = [px;py;pz];
wx = 0;  
wy = 0;
wz = 0;
omega = [wx;wy;wz]; 
vx = 10;
vy = 0;
vz = 0;
v = [vx;vy;vz];
phi = deg2rad(0);
theta = deg2rad(0);
psi = deg2rad(0);
lbd = [phi;theta;psi];
x0 = [p; v; lbd; omega];

a = [0, -vx*sin(theta)+vz*cos(theta), -vx*sin(psi)-vy*cos(psi)
     -vy*sin(phi)-vz*cos(phi), 0, vx*cos(psi)-vy*sin(psi)
      vy*cos(phi)-vz*sin(phi), vx*cos(theta)-vz*sin(theta), 0];

b = [(wy*cos(phi)-wz*sin(phi))*tan(theta), (wy*sin(phi)+wz*cos(phi))/((cos(theta))^2), 0
     -wy*sin(phi)-wz*cos(phi), 0, 0
      (wy*cos(phi)-wz*sin(phi))/cos(theta), ((wy*sin(phi)+wz*cos(phi))*tan(theta))/cos(theta), 0];

c = [0, -g_mars, 0
     g_mars*cos(phi), 0, 0
     -g_mars*sin(phi), 0, 0];


A = [ zeros(3), Euler2R(lbd), a, zeros(3)
     zeros(3), zeros(3), c, skew(v)
     zeros(3), zeros(3), b, Euler2Q(lbd)
     zeros(3), zeros(3), zeros(3), zeros(3)];

B=[zeros(3,4)
   1/ml*f1, 1/ml*f2, 1/ml*f3, 1/ml*f4
   zeros(3,4)
   jl^-1*n1, jl^-1*n2, jl^-1*n3, jl^-1*n4];

C = eye(12);

D = zeros(12,4);

sys = ss(A,B,C,D);
y_L = lsim(sys,u,t,x0);

% Controllability, observability and stability
[M,J] = jordan(A),
if rank(ctrb(A,B)) < size(A,1), disp('The system is not controllable');
else 
    disp('The system is controllable');
end
if rank(obsv(A,C)) < size(A,1), disp('The system is not observable');
else 
    disp('The system is observable');
end

figure(9);
plot(t,y_L(:,1),'k',t,y_L(:,2),'r',t,y_L(:,3),'g',t,y_L(:,4),'y',t,y_L(:,5),'b',t,y_L(:,6),'c');
legend('p_x','p_y','p_z','v_x','v_y','v_z')
grid on;
xlabel('t [s]');

% SVD Frequency Analysis
figure(5);
sigma(sys);
title('Sigma Values (Vx = 10, Vz = Vy = 0)');
grid on;

% Checking stability with LQR controller
Q = blkdiag(1,10,100,1,100,10,0.01*eye(3),0.001*eye(3));
R = 0.1*eye(4);
Klqr = lqr(A,B,Q,R);
lbd_CL_lqr = eig(A-B*Klqr); 
X = A-B*Klqr;
CLlqr = ss(X,B,C,D);
if any(real(lbd_CL_lqr) >= 0), disp('CL system with LQR not stable'); else, disp('CL system with LQR is stable'); end


% Simulation - LQR controller
r_LQR = [0; 500; -21; 0; 0; 0; 0; 0; 0; 0; 0; 0] * ones(1, length(t));
x_LQR = zeros(nx, NSim);
u_LQR = zeros(nu, NSim);
x_LQR(:,1) = [0; 0; -21; 0; 10; 0; 5; 0; 0; 0; 0; 0];
y_LQR = zeros(size(C,1), NSim);

for k = 1:NSim
    % Measurements:
    y_LQR(:,k) = C * x_LQR(:,k);

    % Control:
    u_LQR(:,k) = -Klqr * (x_LQR(:,k) - r_LQR(:,k));

    % Simulation:
    x_dot_LQR = A * x_LQR(:,k) + B * u_LQR(:,k);
    xp_LQR = x_LQR(:,k) + Dt * x_dot_LQR;

    if k < NSim
        x_LQR(:,k+1) = xp_LQR;
    end
end

% Plots - LQR
figure(90);
plot(t, u_LQR);
grid on;
xlabel('t [s]');
ylabel('u(t)');
legend('u_1', 'u_2', 'u_3', 'u_4');
title('Control Inputs with LQR Controller');

figure(91);
plot(t, y_LQR(2,:));
grid on;
xlabel('t [s]');
ylabel('p_y');
legend('p_y');
title('Position p_y with LQR Controller');


% H_inf Controller - Definition

A0 = A;
B1 = zeros(12,12);
B2 = B;
W1 = sqrt(Q);
W2 = sqrt(R);
C1 = [W1];
D11 = zeros(12,12);
D12 = [zeros(8,4); W2];
C2 = -eye(12);
D21 = eye(12);
D22 = zeros(12,4);

B0 = [B1, B2];
C0 = [C1; C2];
D0 = [D11, D12; D21, D22];

P = ss(A0, B0, C0, D0);

nmeas = 12;
ncont = 4; 

% tests on P:
if (12 - rank(ctrb(A0,B2))) > 0, disp('A1.1 on P: system uncontrolable'); else disp('A1.1 on P: OK'); end
if (12 - rank(obsv(A0,C2))) > 0, disp('A1.2 on P: system unobservable'); else disp('A1.2 on P: OK'); end
if (size(D12,2) - rank(D12)) > 0, disp('A2.1 on P: D12 is column rank deficient'); else disp('A2.1 on P: OK'); end
if (size(D21,1) - rank(D21)) > 0, disp('A2.2 on P: D21 is row rank deficient'); else disp('A2.1 on P: OK'); end
syms w real; 
Aux1 = [A0 - j*w*eye(size(A0)) , B2 ; C1 , D12];
if (size(Aux1,2) - rank(Aux1)) > 0,  disp('A3 on P: matrix is column rank deficient'); else disp('A3 on P: OK'); end
Aux2 = [A0 - j*w*eye(size(A0)) , B1 ; C2 , D21];
if (size(Aux2,1) - rank(Aux2)) > 0,  disp('A4 on P: matrix is column rank deficient'); else disp('A4 on P: OK'); end

[Kinf,CLinf,gammainf,info_inf] = hinfsyn(P,nmeas,ncont);
poles_CLinf = pole(CLinf);
if any(real(poles_CLinf) >= 0), disp('CL system with Hinf controller not stable'); else, disp('CL system with Hinf controller is stable'); end

% Hinf controller
Dt = 0.01;
t = 0:Dt:60;
r = [0;500;-21;0;0;0;0;0;0;0;0;0]*(t>=0);
NSim = length(t);
nx = 12;
nu = 4;
x = zeros(nx,NSim);
xu = zeros(nx,NSim);
u = zeros(nu,NSim);
x(:,1) = [0;0;-21;0;10;0;5;0;0;0;0;0];


r = [0;500;-21;0;0;0;0;0;0;0;0;0] * ones(1, length(t));
y1 = zeros(size(C, 1), length(t));

for k = 1:NSim

    % Measurements:
    y1(:,k) = C*x(:,k);

    % Control:
    
    v = -[(r(:,k)-y1(:,k))];
    u(:,k) = Kinf.C*v; 
    
    
    % Simulation:
    x_dot = A*x(:,k) + B*u(:,k); 
    xp = x(:,k) + Dt*x_dot; 
    if k < NSim
        x(:,k+1) = xp;
    end
end

figure(34);
plot(t,u);
grid on;
xlabel('t [s]');
ylabel('u(t)');
legend('u_1', 'u_2', 'u_3', 'u_4');
title('Control Inputs with Hinf Controller');

figure(35);
plot(t,y1(2,:)); 
xlabel('t [s]');
ylabel('p_y');
legend('p_y');
title('Position p_y with Hinf Controller');
grid on;

figure(113);
sigmaplot(CLlqr,'k',Kinf);
grid on;
title('Sigma Values - CL (modelo 4)');
legend('LQR','H\infty');

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