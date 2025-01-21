clear
clc
close all

global J1 J2 m1 m2 r d l k b g sigma_0 sigma_1 sigma_2 theta_s_dot Ts Tc 
global thetaDesired1 thetaDesired2 thetaDesired1_dot thetaDesired2_dot thetaDesired1_doubledot thetaDesired2_doubledot 
global Phi_g1 Phi_f1 Phi_g2 Phi_f2
global u1_values u2_values time_values

J1 = 0.5; J2 = 0.625; m1 = 2; m2 = 2.5;
r = 0.5; d = 0.5; l = 0.5; k = 150; b = 1; g = 9.81;
sigma_0 = 1; sigma_1 = 1; sigma_2 = 1; theta_s_dot = 0.1;
Ts = 2; Tc = 1;

u1_values = [];
u2_values = [];
time_values = [];

%% Define Desired Trajectories
thetaDesired1 = @(t) pi * sin(2 * pi * t) * (1/6);
thetaDesired2 = @(t) pi * sin(pi * t) * (1/4);
thetaDesired1_dot = @(t) (pi^2/3) * cos(2 * pi * t);
thetaDesired2_dot = @(t) (pi^2/4) * cos(pi * t);
thetaDesired1_doubledot = @(t) -2/3 * pi^3 * sin(2 * pi * t);
thetaDesired2_doubledot = @(t) -(pi^3/4) * sin(pi * t);

% thetaDesired1 = @(t) (pi/6) * cos(2 * pi * t) + (pi/4) * cos(3 * pi * t);
% thetaDesired1_dot = @(t) - (pi^2/3) * sin(2 * pi * t) - (3 * pi^2 / 4) * sin(3 * pi * t);
% thetaDesired1_doubledot = @(t) - (pi^3/3) * cos(2 * pi * t) - (9 * pi^3 / 4) * cos(3 * pi * t);
% 
% thetaDesired2 = @(t) (pi/4) * cos(pi * t) + (pi/2) * cos(5 * pi * t);
% thetaDesired2_dot = @(t) - (pi^2/4) * sin(pi * t) - (5 * pi^2 / 2) * sin(5 * pi * t);
% thetaDesired2_doubledot = @(t) - (pi^3/4) * cos(pi * t) - (25 * pi^3 / 2) * cos(5 * pi * t);

% thetaDesired1 = @(t) (pi/6) * cos(2 * pi * t) + (pi/4) * cos(3 * pi * t);
% thetaDesired1_dot = @(t) - (pi^2/3) * sin(2 * pi * t) - (3 * pi^2 / 4) * sin(3 * pi * t);
% thetaDesired1_doubledot = @(t) - (pi^3/3) * cos(2 * pi * t) - (9 * pi^3 / 4) * cos(3 * pi * t);
% 
% thetaDesired2 = @(t) (pi/4) * cos(pi * t) + (pi/2) * cos(5 * pi * t);
% thetaDesired2_dot = @(t) - (pi^2/4) * sin(pi * t) - (5 * pi^2 / 2) * sin(5 * pi * t);
% thetaDesired2_doubledot = @(t) - (pi^3/4) * cos(pi * t) - (25 * pi^3 / 2) * cos(5 * pi * t);


%% Define Basis Functions 
% Define the centers as 2D coordinates

max_theta = pi/2; % this is the maximum expected value for theta
max_thetadot = 10;% this is the maximum expected value for thetadot

centers_f1 = [-max_theta, -max_thetadot; -max_theta, max_thetadot; max_theta, -max_thetadot; max_theta, max_thetadot];
distances = zeros(4, 4);
for i = 1:4
    for j = i+1:4
        distances(i,j) = sqrt((centers_f1(i,1) - centers_f1(j,1))^2 + (centers_f1(i,2) - centers_f1(j,2))^2);
    end
end

mean_distance = mean(distances(:));
sigma_f1 = mean_distance / sqrt(2);

centers_g1 = [-1, 0; 0, 1; 1, -1; -0.5, 0.5];
sigma_g1 = 0.5;

centers_f2 = centers_f1;
sigma_f2 = sigma_f1;

centers_g2 = centers_g1;
sigma_g2 = 0.75;

% Define the Phi functions as functions of x_ and y_
Phi_g1 = @(x_, y_) exp(-sum((repmat([x_, y_], size(centers_g1, 1), 1) - centers_g1).^2, 2) / (2 * sigma_g1^2))';
Phi_f1 = @(x_, y_) exp(-sum((repmat([x_, y_], size(centers_f1, 1), 1) - centers_f1).^2, 2) / (2 * sigma_f1^2))';
Phi_g2 = @(x_, y_) exp(-sum((repmat([x_, y_], size(centers_g2, 1), 1) - centers_g2).^2, 2) / (2 * sigma_g2^2))';
Phi_f2 = @(x_, y_) exp(-sum((repmat([x_, y_], size(centers_f2, 1), 1) - centers_f2).^2, 2) / (2 * sigma_f2^2))';

%% Define TimeSpan for Integration
tmax = 2;
tspan = [0 tmax];

%% Define Initial Values and Run ODE
xInitial = zeros(1, 22);
[t, x] = ode45(@odefun, tspan, xInitial);

%% 
% Evaluate the functions
theta1_values = thetaDesired1(t);
theta2_values = thetaDesired2(t);

e1 = x(:, 1) - theta1_values;
e2 = x(:, 3) - theta2_values;

%% Plots
plot(t, x(:, 1), 'b', 'LineWidth', 2);
hold on
plot(t, theta1_values, '--r', 'LineWidth', 2);
legend('Real', 'Desired', 'Interpreter', 'latex');
title('Tracking for System 1', 'Interpreter', 'latex');
%print('plot1_1', '-depsc');
figure

plot(t, x(:, 2), 'b', 'LineWidth', 2);
title('$\dot{\theta_1}$', 'Interpreter', 'latex');
%print('plot2_1', '-depsc');
figure

plot(t, x(:, 3), 'b', 'LineWidth', 2);
hold on
plot(t, theta2_values, '--r', 'LineWidth', 2);
legend('Real', 'Desired', 'Interpreter', 'latex');
title('Tracking for System 2', 'Interpreter', 'latex');
%print('plot3_1', '-depsc');
figure

plot(t, x(:, 4), 'b', 'LineWidth', 2);
title('$\dot{\theta_2}$', 'Interpreter', 'latex');
%print('plot4_1', '-depsc');
figure

plot(t, e1, 'b', 'LineWidth', 2);
title('Tracking Error for Dynamic System 1: $\theta_1 - \theta_{d1}$', 'Interpreter', 'latex');
%print('plot5_1', '-depsc');
figure

plot(t, e2, 'b', 'LineWidth', 2);
title('Tracking Error for Dynamic System 2: $\theta_2 - \theta_{d2}$', 'Interpreter', 'latex');
%print('plot6_1', '-depsc');
figure

plot(t, x(:, 5), 'b', 'LineWidth', 2);
hold on
plot(t, x(:, 6), 'r', 'LineWidth', 2);
legend('$\tau_1$', '$\tau_2$', 'Interpreter', 'latex'); 
%print('plot7_1', '-depsc');
figure

plot(t, x(:, 7), 'b', 'LineWidth', 2)
hold on
plot(t, x(:, 8), 'r', 'LineWidth', 2)
hold on
plot(t, x(:, 9), 'g', 'LineWidth', 2)
hold on
plot(t, x(:, 10), 'm', 'LineWidth', 2)
legend('$\theta_{f1}$', '$\theta_{f2}$', '$\theta_{f3}$', '$\theta_{f4}$', 'Interpreter', 'latex');
title('Weights Estimations for $f_1(x)$', 'Interpreter','latex')
%print('plot8_1', '-depsc');
figure

plot(t, x(:, 11), 'b', 'LineWidth', 2)
hold on
plot(t, x(:, 12), 'r', 'LineWidth', 2)
hold on
plot(t, x(:, 13), 'g', 'LineWidth', 2)
hold on
plot(t, x(:, 14), 'm', 'LineWidth', 2)
legend('$\theta_{g1}$', '$\theta_{g2}$', '$\theta_{g3}$', '$\theta_{g4}$', 'Interpreter', 'latex');
title('Weights Estimations for $g_1(x)$', 'Interpreter','latex')
%print('plot9_1', '-depsc');
figure

plot(t, x(:, 15), 'b', 'LineWidth', 2)
hold on 
plot(t, x(:, 16), 'r', 'LineWidth', 2)
hold on
plot(t, x(:, 17), 'g', 'LineWidth', 2)
hold on
plot(t, x(:, 18), 'm', 'LineWidth', 2)
legend('$\theta_{f1}$', '$\theta_{f2}$', '$\theta_{f3}$', '$\theta_{f4}$', 'Interpreter', 'latex');
title('Weights Estimations for $f_2(x)$', 'Interpreter','latex')
%print('plot10_1', '-depsc');
figure

plot(t, x(:, 19), 'b', 'LineWidth', 2)
hold on
plot(t, x(:, 20), 'r', 'LineWidth', 2)
hold on
plot(t, x(:, 21), 'g', 'LineWidth', 2)
hold on
plot(t, x(:, 22), 'm', 'LineWidth', 2)
legend('$\theta_{g1}$', '$\theta_{g2}$', '$\theta_{g3}$', '$\theta_{g4}$', 'Interpreter', 'latex');
title('Weights Estimations for $g_2(x)$', 'Interpreter','latex');
%print('plot11_1', '-depsc');
figure;

plot(time_values, u1_values, 'b', 'LineWidth', 2);
hold on;
plot(time_values, u2_values, 'r', 'LineWidth', 2);
legend('$u_1$', '$u_2$', 'Interpreter', 'latex');
title('Control Inputs $u_1$ and $u_2$', 'Interpreter', 'latex');
xlabel('Time (s)');
ylabel('Control Inputs');
%print('plot12_1', '-depsc');

