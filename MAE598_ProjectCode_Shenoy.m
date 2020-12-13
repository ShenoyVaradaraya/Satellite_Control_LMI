clc;clear;close all;

Ix = 1030.17;
Iy = 3015.65;
Iz = 3030.43;
w0 = 7.292115e-5;
L1 = eye(3);
L2 = eye(3);
A_e41 = -(4*(w0^2)*(Iy-Iz)/Ix);
A_e46 = -w0*(Iy-Iz-Ix)/Ix;
A_e52 = -(3*(w0^2)*(Ix-Iz)/Iy);
A_e63 = -((w0^2)*(Iy-Ix)/Iz);
A_e64 = w0*(Iy-Ix-Iz)/Iz;
A = [ 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0  1;A_e41 0 0 0 0 A_e46;0 A_e52 0 0 0 0; 0 0 A_e63 A_e64 0 0];
B1 = [ zeros(3);inv(Ix) 0 0; 0 inv(Iy) 0; 0 0 inv(Iz)];
B2 = B1;
C1 = 10e-3 * [ -4*(w0^2)*(Iy-Iz) 0 0 0 0 -w0*(Iy-Ix-Iz); 0 -(3*(w0^2)*(Ix-Iz)) 0 0 0 0; 0 0 -(w0^2)*(Iy-Ix) w0*(Iy-Ix-Iz) 0 0];
C2 = [ eye(3) zeros(3)];
D1 = 1e-3*L1;
D2 = 1e-3*L2;
%% H-infinity LMI (mosek)
[~,n] = size(A);
[~,b1_2] = size(B1);
[c1,~] = size(C1);
X = sdpvar(n);
W = sdpvar(c1,n);
gamma = sdpvar(1);
infinity_LMI = [[((A*X+B1*W)' + A*X + B1*W) B2 (C1*X+D1*W)'; B2' -gamma*eye(b1_2) D2' ; (C1*X + D1*W) D2 -gamma*eye(c1)] <= 0];
infinity_LMI = [infinity_LMI,X>=0];
opt = sdpsettings('solver','mosek','verbose',0);
optimize(infinity_LMI,gamma,opt);
opt_gamma = value(gamma);
K = value(W)*inv(value(X));
%% Closed Loop system 
tinit = 0;
tinc  = 0.005;
tfin  = 20.0;
t     = [tinit:tinc:tfin]';   
u0    = 0*[t t t];
u     = [sin(t/200) cos(t/200) sin(t/250)];  % Vector of uniformly spaced time points
x0 = [0.05 0.05 0.05 0.01 0.01 0.01]; 
ty0 = lsim(ss((A+B1*K),B2,(C1+D1*K),D2), u0,t,x0);
ty1 = lsim(ss((A+B1*K),B2,(C1+D1*K),D2), u,t,x0);
figure()
subplot(3,1,1)
hold on
plot(t,ty0(:,1))
plot(t,ty1(:,1))
ylabel('Output z_{\infty 1}');
xlabel('Time(sec)');
legend('d = 0','d = [sin(t/200) cos(t/200) sin(t/250)]^{T}')
title('H_{\infty} Control');
grid
hold off
subplot(3,1,2)
hold on
plot(t,ty0(:,2))
plot(t,ty1(:,2))
ylabel('Output z_{\infty 2}');
xlabel('Time(sec)');
legend('d = 0','d = [sin(t/200) cos(t/200) sin(t/250)]^{T}')
grid
hold off
subplot(3,1,3)
hold on
plot(t,ty0(:,3))
plot(t,ty1(:,3))
ylabel('Output z_{\infty 3}');
xlabel('Time(sec)');
legend('d = 0','d = [sin(t/200) cos(t/200) sin(t/250)]^{T}')
grid 
hold off
pause
%% H2 LMI(mosek)
[~,n] = size(A);
[c1,~] = size(C2);
X2 = sdpvar(n);
W2 = sdpvar(c1,n);
Z = sdpvar(c1,c1);
rho = 0.003;
H2_LMI = [(A*X2 + B1*W2 + (A*X2 + B1*W2)' + B2*B2') <= 0];
H2_LMI = [H2_LMI,[-Z (C2*X2); X2*C2' -X] <= 0];
H2_LMI = [H2_LMI,trace(Z) <= rho];
H2_LMI = [H2_LMI,X2>= 0];
optimize(H2_LMI,rho,opt);
opt_gamma2 = sqrt(value(rho));
K2 = value(W2)*inv(value(X2));
controller2 = ss(K2);
%% Closed Loop System
tinit = 0;
tinc  = 0.005;
tfin  = 20.0;
t     = [tinit:tinc:tfin]';    
u0    = 0*[t t t];
u2     = 10*[sin(t) cos(t) sin(2*t)];  % Vector of uniformly spaced time points
x0 = [0.05 0.05 0.05 0.01 0.01 0.01]; 
ty0 = lsim(ss((A+B1*K2),B2,C2,0), u0,t,x0);
ty2 = lsim(ss((A+B1*K2),B2,C2,0), u2,t,x0);
figure()
subplot(3,1,1)
hold on
plot(t,ty0(:,1))
plot(t,ty2(:,1))
ylabel('Output z_{21}');
xlabel('Time(sec)');
legend('d = 0','d = [sin(t) cos(t) sin(2t)]^{T}')
title('H_{2} Control');
grid
hold off
subplot(3,1,2)
hold on
plot(t,ty0(:,2))
plot(t,ty2(:,2))
ylabel('Output z_{22}');
xlabel('Time(sec)');
legend('d = 0','d = [sin(t) cos(t) sin(2t)]^{T}')
grid
hold off
subplot(3,1,3)
hold on
plot(t,ty0(:,3))
plot(t,ty2(:,3))
ylabel('Output z_{23}');
xlabel('Time(sec)');
legend('d = 0','d = [sin(t) cos(t) sin(2t)]^{T}')
grid 
hold off
pause
%% mixed H2/Hinf 
L = [ 2 0; 0 -6];
M = [ 1 0; 0 -1];
X = sdpvar(n);
W = sdpvar(c1,n);
Z = sdpvar(c1,c1);
gamma = sdpvar(1);
rho = sdpvar(1);
h2hinf_LMI = [ [ -Z C2*X; X*C2' -X] <= 0];
h2hinf_LMI = [h2hinf_LMI, trace(Z) <=rho];
h2hinf_LMI = [h2hinf_LMI, A*X+B1*W+(A*X+B1*W)'+B2*B2'<=0];
h2hinf_LMI = [h2hinf_LMI, kron(L,X)+kron(M,(A*X+B1*W))+kron(M',(A*X+B1*W)') <=0];
h2hinf_LMI =[h2hinf_LMI,[((A*X+B1*W)' + A*X + B1*W) B1 (C1*X+D1*W)'; 
                    B1' -gamma*eye(b1_2) D2' ; 
                    (C1*X + D1*W) D2 -gamma*eye(c1)] <= 0];           
optimize(h2hinf_LMI,gamma+rho,opt);
K2inf=value(W)*inv(value(X));
%% Closed Loop System 
tinit = 0;
tinc  = 0.005;
tfin  = 20;
t     = [tinit:tinc:tfin]';    
u3     = 1e-3*[sin(w0*t) sin(w0*t) sin(w0*t)];  % Vector of uniformly spaced time points
x0 = [0.05 0.05 0.05 0.01 0.01 0.01]; 
ty1a = lsim(ss((A+B1*K),B2,(C2+D1*K),D2), u3,t,x0);
ty3 = lsim(ss((A+B1*K2inf),B2,(C2+D1*K2inf),D2), u3,t,x0);
figure()
subplot(3,1,1)
hold on
plot(t,ty1a(:,1))
plot(t,ty3(:,1))
ylabel('Component 1');
xlabel('Time(sec)');
legend('Ordinary Feedback Control','Mixed H_{2}/H_{\infty} Feedback Control')
title('Mixed H_{2}/H_{\infty} Control');
grid
hold off
subplot(3,1,2)
hold on
plot(t,ty1a(:,2))
plot(t,ty3(:,2))
ylabel('Component 2');
xlabel('Time(sec)');
legend('Ordinary Feedback Control','Mixed H_{2}/H_{\infty} Feedback Control')
grid
hold off
subplot(3,1,3)
hold on
plot(t,ty1a(:,3))
plot(t,ty3(:,3))
ylabel('Component 3');
xlabel('Time(sec)');
legend('Ordinary Feedback Control','Mixed H_{2}/H_{\infty} Feedback Control')
grid
hold off
pause
ty2a = lsim(ss((A+B1*K2),B2,C2,0), u3,t,x0);
ty4 = lsim(ss((A+B1*K2inf),B2,C2,0), u3,t,x0);
figure()
subplot(3,1,1)
hold on
plot(t,ty2a(:,1))
plot(t,ty4(:,1))
ylabel('Component 1');
xlabel('Time(sec)');
legend('Ordinary Feedback Control','Mixed H_{2}/H_{\infty} Feedback Control')
title('Mixed H_{2}/H_{\infty} Control');
grid
hold off
subplot(3,1,2)
hold on
plot(t,ty2a(:,2))
plot(t,ty4(:,2))
ylabel('Component 2');
xlabel('Time(sec)');
legend('Ordinary Feedback Control','Mixed H_{2}/H_{\infty} Feedback Control')

grid
hold off
subplot(3,1,3)
hold on
plot(t,ty2a(:,3))
plot(t,ty4(:,3))
ylabel('Component 3');
xlabel('Time(sec)');
legend('Ordinary Feedback Control','Mixed H_{2}/H_{\infty} Feedback Control')
grid 
hold off
pause
%% H-infinity state feedback synthesis for interval uncertainty
clear Ix Iy Iz A B1 B2 C1 A_e41 A_e46 A_e52 A_e63 A_e64
Ix = ureal('Ix',1030.17,'Range',[1030.17 1960.21]);
Iy = ureal('Iy',3015.65,'Range',[3015.65 3260.59 ]);
Iz = ureal('Iz',3030.43,'Range',[3030.43 3385.48]);
A_e41 = -(4*(w0^2)*(Iy-Iz)/Ix);
A_e46 = -w0*(Iy-Iz-Ix)/Ix;
A_e52 = -(3*(w0^2)*(Ix-Iz)/Iy);
A_e63 = -((w0^2)*(Iy-Ix)/Iz);
A_e64 = w0*(Iy-Ix-Iz)/Iz;
A = [ 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0  1;A_e41 0 0 0 0 A_e46;0 A_e52 0 0 0 0; 0 0 A_e63 A_e64 0 0];
B1 = [ zeros(3);inv(Ix) 0 0; 0 inv(Iy) 0; 0 0 inv(Iz)];
B2 = B1;
C1 = 10e-3 * [ -4*(w0^2)*(Iy-Iz) 0 0 0 0 -w0*(Iy-Ix-Iz); 0 -(3*(w0^2)*(Ix-Iz)) 0 0 0 0; 0 0 -(w0^2)*(Iy-Ix) w0*(Iy-Ix-Iz) 0 0];
P = ss(A,[B1 B2],C1,[D1 D2]);
[Kint,Clint,intgam] = hinfsyn(P,2,2);
impulse(Clint)
%--------------------------------------------------------------------------%
