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
zeta_chi = 1;


%tuning variables
kp_fi = -2;
ki_fi = -0;
kd_fi = (2*zeta_fi*omega_fi-alpha_fi1)/alpha_fi2;
kp_chi = 2*zeta_chi*omega_chi*V_g/g;
ki_chi = omega_chi^2*V_g/g;

fi_old = fi;
chi_old = chi;
delta_a_old = delta_a;

sim('course_hold_AP_2015');

%% Plot

figure(1)
subplot(3,1,1)
plot(fi.time,fi.data*180/pi,'b')
hold on
plot(fi_old.time,fi_old.data*180/pi,'r')
legend('fi')
subplot(3,1,2)
plot(chi.time,chi.data*180/pi,'b')
hold on
plot(chi_old.time,chi_old.data*180/pi,'r')
hold on
plot(chi_ref.time,chi_ref.data,'g')
legend('chi','chi prev','chi ref')
subplot(3,1,3)
plot(delta_a.time,delta_a.data*180/pi,'b')
hold on
plot(delta_a_old.time,delta_a_old.data*180/pi,'r')
legend('delta_a')






