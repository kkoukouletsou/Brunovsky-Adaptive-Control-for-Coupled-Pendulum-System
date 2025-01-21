clear
clc
close all

x_f = linspace(-20, 20, 25);
y_f = linspace(-20, 20, 25);
[X_f, Y_f] = meshgrid(x_f, y_f);

x_g = linspace(-2, 2, 25);
y_g = linspace(-4, 4, 25);
[X_g, Y_g] = meshgrid(x_g, y_g);

max_theta = pi;
max_thetadot = 10;

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

Phi_f1_vals = zeros(size(X_f, 1), size(X_f, 2), size(centers_f1, 1));
Phi_g1_vals = zeros(size(X_g, 1), size(X_g, 2), size(centers_g1, 1));
Phi_f2_vals = zeros(size(X_f, 1), size(X_f, 2), size(centers_f2, 1));
Phi_g2_vals = zeros(size(X_g, 1), size(X_g, 2), size(centers_g2, 1));

for i = 1:size(centers_f1, 1)
    Phi_f1_vals(:, :, i) = exp(-((X_f - centers_f1(i, 1)).^2 + (Y_f - centers_f1(i, 2)).^2) / (2 * sigma_f1^2));
end

for i = 1:size(centers_g1, 1)
    Phi_g1_vals(:, :, i) = exp(-((X_g - centers_g1(i, 1)).^2 + (Y_g - centers_g1(i, 2)).^2) / (2 * sigma_g1^2));
end

for i = 1:size(centers_f2, 1)
    Phi_f2_vals(:, :, i) = exp(-((X_f - centers_f2(i, 1)).^2 + (Y_f - centers_f2(i, 2)).^2) / (2 * sigma_f2^2));
end

for i = 1:size(centers_g2, 1)
    Phi_g2_vals(:, :, i) = exp(-((X_g - centers_g2(i, 1)).^2 + (Y_g - centers_g2(i, 2)).^2) / (2 * sigma_g2^2));
end

figure;
for i = 1:size(centers_f1, 1)
    subplot(2, 2, i);
    surf(X_f, Y_f, Phi_f1_vals(:, :, i));
    title(['$\Phi_{f1}$ Center ', num2str(i)], 'Interpreter', 'latex');
    xlabel('x');
    ylabel('y');
    zlabel('$\Phi$', 'Interpreter', 'latex');
end

figure;
for i = 1:size(centers_g1, 1)
    subplot(2, 2, i);
    surf(X_g, Y_g, Phi_g1_vals(:, :, i));
    title(['$\Phi_{g1}$ Center ', num2str(i)], 'Interpreter', 'latex');
    xlabel('x');
    ylabel('y');
    zlabel('$\Phi$', 'Interpreter', 'latex');
end

figure;
for i = 1:size(centers_f2, 1)
    subplot(2, 2, i);
    surf(X_f, Y_f, Phi_f2_vals(:, :, i));
    title(['$\Phi_{f2}$ Center ', num2str(i)], 'Interpreter', 'latex');
    xlabel('x');
    ylabel('y');
    zlabel('$\Phi$', 'Interpreter', 'latex');
end

figure;
for i = 1:size(centers_g2, 1)
    subplot(2, 2, i);
    surf(X_g, Y_g, Phi_g2_vals(:, :, i));
    title(['$\Phi_{g2}$ Center ', num2str(i)], 'Interpreter', 'latex');
    xlabel('x');
    ylabel('y');
    zlabel('$\Phi$', 'Interpreter', 'latex');
end
