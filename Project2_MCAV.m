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

%% MODEL 1 & 2

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

C = [eye(3), zeros(3), zeros(3), zeros(3)
     zeros(1,3) zI.' zeros(1,3) zeros(1,3)];

D = zeros(4);

sys = ss(A,B,C,D);
y_L = lsim(sys,u,t,x0);
[V,D,W] = eig(A); 

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
ctrb_modes = W'*B,
obsv_modes = C*V,



%% MODEL 3 & 4

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
   1/m*f1, 1/m*f2, 1/m*f3, 1/m*f4
   zeros(3,4)
   jl^-1*n1, jl^-1*n2, jl^-1*n3, jl^-1*n4];

C = [eye(3), zeros(3), zeros(3), zeros(3)
     zeros(3), eye(3), zeros(3), zeros(3)];

D=zeros(6,4);

sys = ss(A,B,C,D);
y_L = lsim(sys,u,t,x0);


[V,D,W] = eig(A); 

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
ctrb_modes = W'*B,
obsv_modes = C*V,


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