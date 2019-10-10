a_fi1 = 2.87;
a_fi2 = -0.65; % kan være +7.5
delta_max = 30;
epsilon_max = 15;
zeta_fi = 0.707;

omega_fi = sqrt(abs(a_fi2)*delta_max/epsilon_max);

K_pfi = -2;
K_dfi = (2*zeta_fi*omega_fi-a_fi1)/a_fi2;
%K_ifi = 0.6;

h = tf([a_fi2],[1 (a_fi1+a_fi2*K_dfi) (a_fi2*K_pfi) 0]);
g = -h;


figure(1)
rlocus(h,g)
legend('H','-H')
ylim([-1.5 1.5])
xlim([-1.5 1])