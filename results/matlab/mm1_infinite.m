T = readtable('results/mm1_timeseries.csv');

% Graph 1: Average number in system over time
figure;
plot(T.MC, T.Ns_running, 'LineWidth', 1.3);
hold on;
yline(2.571429, 'r--', 'N_s theory = 2.57');
xlabel('Time');
ylabel('N_s');
title('Average Number in System');
grid on;
legend('Simulation', 'Theory');
% Graph 2: Average number in queue over time

figure;
plot(T.MC, T.Nq_running, 'LineWidth', 1.3);
hold on;
yline(1.851429, 'r--', 'N_q theory = 1.85');
xlabel('Time');
ylabel('N_q');
title('Average Number in Queue');
grid on;
legend('Simulation', 'Theory');


% Graph 3: Average system delay over completed customers
figure;
plot(T.completed, T.Ds_running, 'LineWidth', 1.3);
hold on;
yline(14.285714, 'r--', 'D_s theory = 14.29');
xlabel('Completed Customers');
ylabel('D_s');
title('Average System Delay vs. Completions');
grid on;
legend('Simulation', 'Theory');



% Graph 4: Average queue delay over completed customers
figure;
plot(T.completed, T.Dq_running, 'LineWidth', 1.3);
hold on;
yline(10.2857, 'r--', 'D_q theory = 10.29');
xlabel('Completed Customers');
ylabel('D_q');
title('Average Queue Delay vs. Completions');
grid on;
legend('Simulation', 'Theory');
 