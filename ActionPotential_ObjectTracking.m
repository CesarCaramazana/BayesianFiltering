clear all;
clc;
% FitzHugh-Nagumo model



% 2.2.1 Simulation of the model

F = 0.5;
eps = 1/12;
A = 0.7;
B = 0.8;
h = 0.1; %Timestep
sigma_z = h/2;
sigma_s = 4;

u0 = rand(); %U(0,1)
v0 = rand();

u = u0;
v = v0;
m = 2e4;

u_n = zeros(m, 1); %To save component U (state)
v_n = zeros(m, 1); %To save component V (state)
y_n = zeros(m, 1); %To save observations 

x_hat = zeros(2, m); % X estimate [v, u]


N = 100; %Number of particles
u_i = rand(1,N);
v_i = rand(1,N);

for t=(1:m)
    z = sqrt(sigma_z) * randn(); % Noise in state
    s = sqrt(sigma_s) * randn(); % Noise in observation
    
    y = v^2 + s; %Get an observation
    
    %Update
    v = v + h*(v - (v^3)/3 - u + F); %State first component
    u = u + h*eps*(v + A - B*u) + z; %State second component

    %Save values
    v_n(t) = v;
    u_n(t) = u;
    y_n(t) = y;

    %-------------- PARTICLE FILTER -----------

    Z_i = sqrt(sigma_z) * randn(1,N);

    
    u_i_n1 = u_i; % u_i(n-1)
    u_i = u_i + h*eps*(v_i + A - B*u_i) + Z_i;
    v_i = v_i + h * (v_i - (v_i.^3)/3 - u_i_n1 + F);

    x_i = [v_i; u_i];
    
    % Compute weights
    lambda_i = (-1/2*sigma_s) * (y - v_i.^2).^2; %?

    L_n = max(lambda_i);

    w_i = exp(lambda_i - L_n);
    w_i = w_i / sum(w_i);

    % Resampling

    idx = randsample(1:N, N, true, w_i);

    X_old = x_i;
    x_i(:, 1:N) = X_old(:, idx);

    v_i = x_i(1,:);
    u_i = x_i(2,:);

    % d) Estimate the state

    x_h_Particle = sum(w_i.*x_i, 2);

    % Save the estimate

    x_hat(:,t) = x_h_Particle; %Estimate


end

%%

%Plots

figure,
subplot(2,1,1)
plot(v_n, 'k-', 'LineWidth',3.0); 
hold on; 
%plot(y_n, 'r--');
plot(x_hat(1,:), 'b--');
title("V: Neuron membrane potential");
legend("Real", "Estimate");
xlabel("Time");
ylabel("Volts");

subplot(2,1,2)
plot(u_n, 'k-', 'LineWidth',1.5);
hold on;
plot(x_hat(2, :), 'b--');
title("U: Recovery current");
xlabel("Time");
ylabel("Milliamps");
legend("Real", "Estimate");

figure,
plot(v_n, u_n, 'r-');
hold on;
plot(x_hat(1,:), x_hat(2,:), 'b--');
xlabel("Volts");
ylabel("Milliamps");
title("Potential vs Current");
legend("Real", "Estimate");
grid on;


figure,
subplot(3,1,1);
plot(v_n(1:1000), 'k-', 'LineWidth',1.0); 
xlabel("Time");
title("Real V");
ylabel("Volts");

subplot(3,1,2);
plot(x_hat(1,1:1000), 'b-');
title("Estimate");
xlabel("Time");
ylabel("Volts");

subplot(3,1,3);
plot(y_n(1:1000), 'r-');
title("Observations");
xlabel("Time");
ylabel("Volts");

%%

% Normalized squared errors

v_hat = x_hat(1,:)';
u_hat = x_hat(2,:)';

err_v = sum((v_n - v_hat).^2) / sum(v_n.^2);
err_u = sum((u_n - u_hat).^2) / sum(u_n.^2);

errors = [err_v, err_u];

axs = categorical({'Error V', 'Error U'});
figure,
bar(axs, errors);
title("Normalized average squared errors");
