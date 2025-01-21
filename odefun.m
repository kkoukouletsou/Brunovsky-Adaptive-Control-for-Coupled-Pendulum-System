function dydt = odefun(t, y)

theta1 = y(1);
theta1_dot = y(2);
theta2 = y(3);
theta2_dot = y(4);
t1 = y(5);
t2 = y(6);
fWeightVector1_hat = y(7:10); 
gWeightVector1_hat = y(11:14);
fWeightVector2_hat = y(15:18);
gWeightVector2_hat = y(19:22);

global J1 J2 m1 m2 r d l k b g sigma_0 sigma_1 sigma_2 theta_s_dot Ts Tc
global thetaDesired1 thetaDesired2 thetaDesired1_dot thetaDesired2_dot thetaDesired1_doubledot thetaDesired2_doubledot 
global Phi_g1 Phi_f1 Phi_g2 Phi_f2
global u1_values u2_values time_values


% Define Parameters
lamda1 = 50; 
lamda2 = 90;
n_a1 = 1; n_a2 = 1; n_b1 = 1; n_b2 = 1; delta1 = 1; delta2 = 1; k1 = 50; k2 = 50; n1 = 1; n2 = 1;
a1 = 1; a2 = 2;
sigma_f1 = 5;
sigma_f2 = 3;
sigma_g1 = 4;
sigma_g2 = 2;
Gammaf1 = 5 * eye(4);
Gammaf2 = 5 * eye(4);
Gammag1 = 3 * eye(4);
Gammag2 = 3 * eye(4);
disp(t)

% Controller 
E1 = theta1_dot - thetaDesired1_dot(t) + lamda1 * (theta1 - thetaDesired1(t));
E2 = theta2_dot - thetaDesired2_dot(t) + lamda2 * (theta2 - thetaDesired2(t));

v1 = - thetaDesired1_doubledot(t) + a1 * (theta1_dot - thetaDesired1_dot(t));
v2 = - thetaDesired2_doubledot(t) + a2 * (theta2_dot - thetaDesired2_dot(t)); 

u_b1 = fWeightVector1_hat.' * Phi_f1(theta1, theta1_dot).' + v1 + k1 * E1 + n1 * E1;
u_a1 = - ((gWeightVector1_hat.' * Phi_g1(theta1, theta1_dot).')/(gWeightVector1_hat.' * Phi_g1(theta1, theta1_dot).' + delta1)) * u_b1; 

u_b2 = fWeightVector2_hat.' * Phi_f2(theta2, theta2_dot).' + v2 + k2 * E2 + n2 * E2;
u_a2 = - ((gWeightVector2_hat.' * Phi_g2(theta2, theta2_dot).')/(gWeightVector2_hat.' * Phi_g2(theta2, theta2_dot).' + delta2)) * u_b2;

u1 = u_a1 - (n_a1 * norm(u_a1)^2 + n_b1 * norm(u_b1)^2) * E1; 
u2 = u_a2 - (n_a2 * norm(u_a2)^2 + n_b2 * norm(u_b2)^2) * E2; 

% COMMENT TO SEE CONTROLLER OUTPUT WITHOUT SATURATION
% saturation_value = 20;
% u1 = max(min(u1, saturation_value), -saturation_value); 
% u2 = max(min(u2, saturation_value), -saturation_value);

u1_values = [u1_values; u1];
u2_values = [u2_values; u2];
time_values = [time_values; t];

fWeightVector1_hat_dot = -sigma_f1 * fWeightVector1_hat + E1 * Gammaf1 * Phi_f1(theta1, theta1_dot).';
gWeightVector1_hat_dot = -sigma_g1 * gWeightVector1_hat + E1 * Gammag1 * Phi_g1(theta1, theta1_dot).';

fWeightVector2_hat_dot = -sigma_f2 * fWeightVector2_hat + E2 * Gammaf2 * Phi_f2(theta2, theta2_dot).';
gWeightVector2_hat_dot = -sigma_g2 * gWeightVector2_hat + E2 * Gammag2 * Phi_g2(theta2, theta2_dot).';

% System Dynamics
num = r/2 * (cos(theta2) - cos(theta1));
den = d + r/2 * (sin(theta1) - sin(theta2));
theta = atan(num/den); 

t1dot = theta1_dot - sigma_0 * (abs(theta1_dot)/(Tc + (Ts - Tc)*exp(-abs(theta1_dot/theta_s_dot)))) * t1;
T1 = sigma_0 * t1 + sigma_1 *  t1dot + sigma_2 * theta1_dot;

t2dot = theta2_dot - sigma_0 * (abs(theta2_dot)/(Tc + (Ts - Tc)*exp(-abs(theta2_dot/theta_s_dot)))) * t2;
T2 = sigma_0 * t2 + sigma_1 *  t2dot + sigma_2 * theta2_dot;

x = sqrt(d^2 + d * r * (sin(theta1) - sin(theta2)) + r^2 * ((1-cos(theta2 - theta1))/(2)));
var = 1/(2 * sqrt(d^2 + d * r * (sin(theta1) - sin(theta2)) + r^2 * (1-cos(theta2 - theta1))/(2)));
x_dot = var * (d * r * (cos(theta1) * theta1_dot - cos(theta2) * theta2_dot) + r^2/2 * sin(theta2 - theta1) * (theta2_dot - theta1_dot));
F = k * (x - l) + b * x_dot;


dydt(1) = theta1_dot;
dydt(2) = (m1 * g * r) * (1/J1) * sin(theta1) - (0.5 * F * r) * cos(theta1 - theta) * (1/J1) - T1/J1 + u1/J1;
dydt(3) = theta2_dot;
dydt(4) = (m2 * g * r) * (1/J2) * sin(theta2) + (0.5 * F * r) * cos(theta2 - theta) * (1/J2) - T2/J2 + u2/J2; 
dydt(5) = t1dot;
dydt(6) = t2dot;
dydt(7:10) = fWeightVector1_hat_dot; 
dydt(11:14) = gWeightVector1_hat_dot; 
dydt(15:18) = fWeightVector2_hat_dot; 
dydt(19:22) = gWeightVector2_hat_dot;

dydt = dydt.';

end