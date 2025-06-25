% 1D Plasticity Solver with Nonlinear Hardening
clear all; close all; clc;

% Material Parameters
E = 210000;          % Young's modulus [MPa]
sigma_y0 = 85;       % Initial yield stress [MPa]
h1 = 10.0;           % Hardening parameter 1
h2 = 300.0;          % Hardening parameter 2

% Loading Parameters
num_steps = 1000;    % Number of time steps
strain_max = 0.001;  % Maximum strain amplitude

% Create cyclic loading pattern starting from 0 to +strain_max, 
% then to -strain_max and back
t = linspace(0, 4*pi, num_steps);
strain = strain_max * sin(t);  % Sinusoidal loading pattern

% Initialize arrays
sigma = zeros(1, num_steps);      % Stress
epsilon_p = zeros(1, num_steps);  % Plastic strain
alpha = zeros(1, num_steps);      % Internal variable

% Time integration loop
for k = 2:num_steps
    % Current total strain
    epsilon_k = strain(k);
    
    % Trial state (elastic predictor)
    sigma_trial = E * (epsilon_k - epsilon_p(k-1));
    
    % Current yield stress
    sigma_y = sigma_y0 + h1 * (1 - exp(-h2 * alpha(k-1)));
    
    % Check yield condition
    f_trial = abs(sigma_trial) - sigma_y;
    
    if f_trial <= 0
        % Elastic step
        sigma(k) = sigma_trial;
        epsilon_p(k) = epsilon_p(k-1);
        alpha(k) = alpha(k-1);
    else
        % Plastic step - Newton-Raphson iteration
        delta_gamma = 0;
        tol = 1e-8;
        max_iter = 100;
        
        for iter = 1:max_iter
            % Current residual
            sigma_current = sigma_trial - E * delta_gamma * sign(sigma_trial);
            sigma_y_current = sigma_y0 + h1 * (1 - exp(-h2 * (alpha(k-1) + delta_gamma)));
            R = abs(sigma_current) - sigma_y_current;
            
            % Check convergence
            if abs(R) < tol
                break;
            end
            
            % Compute tangent and update
            dR = -E - h1 * h2 * exp(-h2 * (alpha(k-1) + delta_gamma));
            delta_gamma = delta_gamma - R/dR;
        end
        
        % Update state variables
        sigma(k) = sigma_trial - E * delta_gamma * sign(sigma_trial);
        epsilon_p(k) = epsilon_p(k-1) + delta_gamma * sign(sigma_trial);
        alpha(k) = alpha(k-1) + delta_gamma;
    end
end

% Plot results
figure('Position', [100 100 800 600])
plot(strain, sigma, 'b-', 'LineWidth', 0.5)
grid on
xlabel('Strain \epsilon', 'FontSize', 12)
ylabel('Stress \sigma [MPa]', 'FontSize', 12)
title('Stress-Strain Response under Cyclic Loading', 'FontSize', 14)
axis tight

% Add some visual improvements
set(gca, 'FontSize', 12)
box on