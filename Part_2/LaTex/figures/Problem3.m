sim_time = 600; % seconds

%constants
g = 9.81;
V_g = 580/3.6;
alpha_fi1 = 2.87;
alpha_fi2 = -0.65;
delta_max = 30;
epsilon_max = 15;
zeta_fi = 0.707;

omega_fi = sqrt(abs(alpha_fi2)*delta_max/epsilon_max);

W_chi = 18; % Between 5 and 10
omega_chi = 1/W_chi*omega_fi;
zeta_chi = 1.8;

%tuning variables
kp_fi = -2;
ki_fi = -0;
kd_fi = (2*zeta_fi*omega_fi-alpha_fi1)/alpha_fi2;
kp_chi = 2*zeta_chi*omega_chi*V_g/g;
ki_chi = omega_chi^2*V_g/g;

compare_to_previous = 1;
% if compare_to_previous
%     fi_old = fi;
%     chi_old = chi;
%     delta_a_old = delta_a;
% end

% Reference signal vector
t = [0 60 140 230 300 380 500 600]';
inputs = [0 10 -5 -10 5 10 20 20]';
chi_ref.time = t;
chi_ref.signals.values = inputs;

A_k = [-0.322, 0.052, 0.028, -1.12;
      0, 0, 1, -0.001;
     -10.6, 0, -2.87, 0.46;
      6.87, 0, -0.04, -0.32];
 
B_k = [0.002; 0; -0.65; -0.02];

C_k = [0 0 1 0;
       0 0 0 1];
 
D_k = [0; 0];
 
sys = ss(A_k, B_k, C_k, D_k);
Ts = 0.01;
sysd = c2d(sys, Ts);


h = 0.01;
Q = h*10^-6*[0.001 0 0 0;
             0     1 0 0;
             0     0 100 0;
             0     0 0   10];
         
R = h*pi/180*[0.2^2 0;
             0   0.2^2];
         
%R = [0.1 0; 0 0.1];
%Making data struct used in Kalman filter in simulink
P_0_apriori = [0.1 0 0 0;
             0 0.001 0 0;
             0 0 0.2*pi/180 0;
             0 0 0 0.2*pi/180];
         
P_0_apriori = Q;
% P_0_apriori = [10000 0 0 0;
%              0 1 0 0;
%              0 0 0.00001 0;
%              0 0 0 100];

E_k = eye(4);
A_d = eye(4) + h*A_k;
B_d = h*B_k;
E_d = E_k;

x_0_apriori = zeros(4,1); %%velg dimensjon
data.A = A_d;
data.B = B_d;
data.C = C_k;
data.Q = Q;
data.R = R;
data.P = P_0_apriori;
data.xhat0 = x_0_apriori;
data.E = E_d; %velg dimensjon her

x_0_apriori = zeros(4,1); %%velg dimensjon
data.A = sysd.A;
data.B = sysd.B;
data.C = sysd.C;
data.Q = Q;
data.R = R;
data.P = P_0_apriori;
data.xhat0 = x_0_apriori;
data.I = eye(4); %velg dimensjon her

prev_run = delta_a;

sim('Problem3_sim.slx');

%% Plot

figure(1); clf;
subplot(2,2,1)
plot(phi.time,phi.data*180/pi,'b')
hold on
plot(phi_est.time,phi_est.data*180/pi,'r')
hold on
legend({'$\phi$','$\phi_{est}$'},'Interpreter','latex','Location','southeast')
title('Roll')
ylabel('Angle [deg]')
set(gca,'FontSize',16)
ylim([-50 50])

subplot(2,2,2)
plot(chi.time,chi.data*180/pi,'b')
hold on
stairs(chi_ref.time,chi_ref.signals.values,'k')
hold on
legend({'$\chi$','$\chi_{ref}$'},'Interpreter','latex','Location','southeast')
title('Course')
ylabel('Angle [deg]')
set(gca,'FontSize',18)
ylim([-15 25])

subplot(2,2,3)
plot(delta_a.time,delta_a.data*180/pi,'b')
hold on
legend({'$\delta_a$'},'Interpreter','latex')
title('Ailerons')
ylabel('Angle [deg]')
xlabel('Time [s]')
set(gca,'FontSize',18)
ylim([-35 35])

subplot(2,2,4)
plot(beta.time,beta.data*180/pi,'b')
hold on
plot(beta_est.time,beta_est.data*180/pi,'r')
hold on
legend({'$\beta$','$\beta_{est}$'},'Interpreter','latex')
title('Sideslip')
ylabel('Angle [deg]')
xlabel('Time [s]')
set(gca,'FontSize',18)
ylim([-0.3 0.3])

% P og R
figure(2); clf;
subplot(2,1,1)
plot(p_meas.time,p_meas.data*180/pi,'g')
hold on
plot(p.time,p.data*180/pi,'b')
hold on
plot(p_est.time,p_est.data*180/pi,'r')
hold on
legend({'$p_{measured}$','$p$','$p_{est}$'},'Interpreter','latex')
title('Roll rate')
ylabel('Angular rate [deg/s]')
set(gca,'FontSize',16)
ylim([-10 10])

subplot(2,1,2)
plot(r_meas.time,r_meas.data*180/pi,'g')
hold on
plot(r.time,r.data*180/pi,'b')
hold on
plot(r_est.time,r_est.data*180/pi,'r')
hold on
legend({'$r_{measured}$','$r$','$r_{est}$'},'Interpreter','latex')
title('Yaw rate')
ylabel('Angular rate [deg/s]')
xlabel('Time [s]')
set(gca,'FontSize',16)
ylim([-3 3])

%% Rudder angles
figure(3); clf;
plot(prev_run.time,prev_run.data*180/pi,'b')
hold on
plot(delta_a.time,delta_a.data*180/pi,'r')
hold on

legend({'$\delta_a, 3 e)$','$\delta_a, 3 f)$'},'Interpreter','latex')
title('Ailerons')
ylabel('Angle [deg]')
xlabel('Time [s]')
set(gca,'FontSize',18)
ylim([-25 25])
