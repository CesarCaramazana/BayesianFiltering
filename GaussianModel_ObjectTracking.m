clear all;
clc;

%===============================================================
%César Caramazana Zarzosa
%Statistical Signal Processing

%November 2022

%===============================================================

% Particle filter

m=100;

T = 0.5; %Timestep
sigma_0 = 0.2; %Noise for prior
sigma_u = 0.1; %Noise in velocity
sigma_w = 0.5; %Noise in position
sigma_l = 16; %Noise in the observations


A = [1, 0, T, 0; 0, 1, 0, T; 0, 0, 1, 0; 0, 0, 0, 1]; %Transition model
Q = [sigma_w, 0, 0, 0; 0, sigma_w, 0, 0; 0, 0, sigma_u, 0; 0, 0, 0, sigma_u]; %Noise in transition model
P0 = [sigma_0, 0, 0, 0;0, sigma_0, 0, 0;0, 0, sigma_0, 0;0, 0, 0, sigma_0]; %Prior
H = [1,-0.3,0,0; -0.2,-1,0,0]; %Observation model
Ln = [sigma_l, 0; 0, sigma_l]; %Noise in observation model

N = 100; %Number of particles
x0 = sqrt(sigma_0) * randn([4, N]); %Prior particles

x_i = x0; %Particles
x = sqrt(P0) * randn([4, 1]); %Initial state

P = P0;
x_h_Kalman = sqrt(P0) * randn([4, 1]); %Initial Kalman estimate
x_h_Particle = sqrt(P0) * randn([4, 1]); %Initial PF estimate

y_n = zeros(2, m); %Matrix to save the observations
x_n = zeros(4, m); %Matrix to save the state
x_hat_PF = zeros(4, m); %Matrix to save the predicted state (Particle filter)
x_hat_KF = zeros(4, m); %Matrix to save the predicted state (Kalman filter)

for t = 1:m

    %Update the state X
    x = A*x + (sqrt(Q)*randn([4, 1])); % x(n) = A*x(n-1) + s(n)

    %Receive an observation
    R_n = sqrt(sigma_l) * randn([2, 1]); %Noise in the observation
    y = H*x + R_n;
    
    % --------- KALMAN FILTER -----------
    % Prediction step
    Pn_ = A*P*A' + Q;
    x_ = A*x_h_Kalman;
   

    % Update step
    S_n = H * Pn_* H' + Ln;
    x_h_Kalman = x_ + Pn_ * H' * (inv(S_n)) * (y - H * x_); %Estimated state
    P = Pn_ - Pn_ * H' * (inv(S_n)) * H * Pn_;


    %----------- PARTICLE FILTER ----------
    
    % a) Sample new particles

    Q_ni = sqrt(Q) * randn([4, N]);

    x_i = A*x_i + Q_ni;

    % b) Compute weights

    lambda_i = (-1/2*sigma_l) * vecnorm(y - H*x_i).^2; 
    L_n = max(lambda_i);

    w_i = exp(lambda_i - L_n); %Weights
    w_i = w_i / sum(w_i); %Normalized

    % c) Resampling

    idx = randsample(1:N, N, true, w_i); %Statistics and Machine Learning Toolbox

    X_old = x_i;
    x_i(:, 1:N) = X_old(:, idx);

    % d) Estimate the state

    x_h_Particle = sum(w_i.*x_i, 2);

    %---------------------------------------------------

    %Save values
    x_n(:, t) = x; %Real state
    y_n(:, t) = y; %Observation
    x_hat_PF(:,t) = x_h_Particle; %Estimate Particle F
    x_hat_KF(:,t) = x_h_Kalman; %Estimate Kalman

    

    
end

%%

%PLOTS
    
% Trajectory (r1, r2)
figure,
plot(x_n(1, 1 : end), x_n(2, 1 : end),'k-', 'LineWidth',2.0);
hold on;
plot(x_hat_PF(1, 1 : end), x_hat_PF(2, 1:end), 'b--', 'LineWidth',1.5);
hold on;
plot(x_hat_KF(1, 1 : end), x_hat_KF(2, 1:end), 'r--', 'LineWidth',1.5);
title("Sample trajectory");
legend("Real", "Particle filter", "Kalman filter");
xlabel("r1");
ylabel("r2");



% Squared errors

R_KF = x_n(1:2, :) - x_hat_KF(1:2, :);
R_PF = x_n(1:2, :) - x_hat_PF(1:2, :);
V_KF = x_n(3:4, :) - x_hat_KF(3:4, :);
V_PF = x_n(3:4, :) - x_hat_PF(3:4, :);

errorR_KF = zeros(2,m);
errorR_PF = zeros(2,m);
errorV_KF = zeros(2,m);
errorV_PF = zeros(2,m);

for n=1:m
    % MSE Kalman filter 
    errorR_KF(1,n) = norm(R_KF(1,n))^2;     
    errorV_KF(1,n) = norm(V_KF(1,n))^2;
    errorR_KF(2,n) = norm(R_KF(2,n))^2;     
    errorV_KF(2,n) = norm(V_KF(2,n))^2;

    % MSE Particle filter 
    errorR_PF(1,n) = norm(R_PF(1,n))^2;     
    errorV_PF(1,n) = norm(V_PF(1,n))^2;
    errorR_PF(2,n) = norm(R_PF(2,n))^2;     
    errorV_PF(2,n) = norm(V_PF(2,n))^2;

end

figure,
subplot(3,2,1);
plot(1:m, x_n(3, 1:m), 'k-');
hold on;
plot(1:m, x_hat_PF(3, 1:m), 'b--');
hold on;
plot(1:m, x_hat_KF(3, 1:m), 'r--');
title("Velocity v1");
legend("Real", "Particle filter", "Kalman filter");
xlabel("Time");

subplot(3,2,2);
plot(1:m, x_n(4, 1:m), 'k-'); 
hold on;
plot(1:m, x_hat_PF(4, 1:m), 'b--');
hold on;
plot(1:m, x_hat_KF(4, 1:m), 'r--');
title("Velocity v2");
legend("Real", "Particle filter", "Kalman filter");
xlabel("Time");



subplot(3,2,3);
plot(errorR_KF(1,:), 'r-');
hold on;
plot(errorR_PF(1,:), 'b-');
title("Squared error r1");
xlabel("Time");
legend("Kalman filter", "Particle filter");

subplot(3,2,4);
plot(errorR_KF(2,:), 'r-');
hold on;
plot(errorR_PF(2,:), 'b-');
title("Squared error  r2");
xlabel("Time");
legend("Kalman filter", "Particle filter");

subplot(3,2,5);
plot(errorV_KF(1,:), 'r-');
hold on;
plot(errorV_PF(1,:), 'b-');
title("Squared error  v1");
xlabel("Time");
legend("Kalman filter", "Particle filter");

subplot(3,2,6);
plot(errorV_KF(2,:), 'r-');
hold on;
plot(errorV_PF(2,:), 'b-');
title("Squared error v2");
xlabel("Time");
legend("Kalman filter", "Particle filter");