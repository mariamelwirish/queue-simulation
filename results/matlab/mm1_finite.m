clear; clc;

T = readtable('mm1_finite_timeseries.csv');


% Parameters
lambda = 0.18;
mu     = 0.25;
K      = 100;          % K waiting, 1 in service => states 0..K+1
rho = lambda / mu;

% Handle rho ~ 1 safely
eps_rho = 1e-12;

if abs(rho - 1) < eps_rho
    % Special case: rho = 1
    P0       = 1/(K + 2);
    Pblock   = P0;                  % P_{K+1}
    Ns_th    = (K + 1)/2;           % mean of 0..K+1
    Nserv_th = 1 - P0;              % prob server busy
    Nq_th    = Ns_th - Nserv_th;
else
    % General case: rho != 1
    P0     = (1 - rho) / (1 - rho^(K + 2));
    Pblock = P0 * rho^(K + 1);

    % Correct numerator with powers K+1 and K+2 (capacity is K+1)
    num = 1 - (K + 2)*rho^(K + 1) + (K + 1)*rho^(K + 2);
    den = (1 - rho) * (1 - rho^(K + 2));      % note single (1 - rho) (after cancel)
    Ns_th = rho * num / den;

    Nserv_th = 1 - P0;
    Nq_th    = Ns_th - Nserv_th;
end

% Use lambda_eff for delays (served throughput)
lambda_eff = lambda * (1 - Pblock);
Ds_th = Ns_th / lambda_eff;
Dq_th = Nq_th / lambda_eff;


% Display theoretical results
fprintf('--- Theoretical (M/M/1/K) ---\n');
fprintf('P_block = %.6f\n', Pblock);
fprintf('N_s = %.6f\n', Ns_th);
fprintf('N_q = %.6f\n', Nq_th);
fprintf('D_s = %.6f\n', Ds_th);
fprintf('D_q = %.6f\n', Dq_th);

% ===============================================
% 1) Average number in system (N_s)
% ===============================================
figure;
plot(T.MC, T.Ns_running, 'b', 'LineWidth', 1.3);
hold on;
yline(Ns_th, 'r--', sprintf('N_s theory = %.2f', Ns_th));
xlabel('Time');
ylabel('N_s');
title('Average Number in System (N_s)');
grid on;
legend('Simulation', 'Theory', 'Location', 'best');

% ===============================================
% 2) Average number in queue (N_q)
% ===============================================
figure;
plot(T.MC, T.Nq_running, 'b', 'LineWidth', 1.3);
hold on;
yline(Nq_th, 'r--', sprintf('N_q theory = %.2f', Nq_th));
xlabel('Time');
ylabel('N_q');
title('Average Number in Queue (N_q)');
grid on;
legend('Simulation', 'Theory', 'Location', 'best');

% ===============================================
% 3) Average system delay (D_s)
% ===============================================
figure;
plot(T.completed, T.Ds_running, 'b', 'LineWidth', 1.3);
hold on;
yline(Ds_th, 'r--', sprintf('D_s theory = %.2f', Ds_th));
xlabel('Completed Customers');
ylabel('D_s');
title('Average System Delay (D_s)');
grid on;
legend('Simulation', 'Theory', 'Location', 'best');

% ===============================================
% 4) Average queue delay (D_q)
% ===============================================
figure;
plot(T.completed, T.Dq_running, 'b', 'LineWidth', 1.3);
hold on;
yline(Dq_th, 'r--', sprintf('D_q theory = %.2f', Dq_th));
xlabel('Completed Customers');
ylabel('D_q');
title('Average Queue Delay (D_q)');
grid on;
legend('Simulation', 'Theory', 'Location', 'best');

% ===============================================
% 5) Blocking probability (P_block)
% ===============================================
figure;
plot(T.MC, T.Pblock_running, 'b', 'LineWidth', 1.3);
hold on;
yline(Pblock, 'r--', sprintf('P_{block} theory = %.6f', Pblock));
xlabel('Time');
ylabel('Blocking Probability');
title('Blocking Probability vs. Time');
grid on;
legend('Simulation', 'Theory', 'Location', 'best');
