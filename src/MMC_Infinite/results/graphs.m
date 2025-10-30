% ---------- User params (must match the sim) ----------
lambda = 0.18;
mu     = 0.25;
c      = 3;

% ---------- Read simulation CSV ----------
T = readtable('mmc_timeseries.csv');  % columns: see header in Java

% Use the logged rows (every K completions) for smooth curves
MC = T.MC;
completed = T.completed;
Ns_run = T.Ns_running;
Nq_run = T.Nq_running;
Ds_run = T.Ds_running;
Dq_run = T.Dq_running;

% Theory from CSV (constant columns) OR compute here
Ns_th_csv = T.Ns_th(end);
Nq_th_csv = T.Nq_th(end);
Ds_th_csv = T.Ds_th(end);
Dq_th_csv = T.Dq_th(end);

% (Optional) Compute theory in MATLAB (Erlang-C)
a   = lambda / mu;
rho = lambda / (c*mu);
if rho >= 1
    error('Unstable: rho >= 1');
end
P0 = 0;
for n=0:c-1
    P0 = P0 + a^n / factorial(n);
end
P0 = 1 / (P0 + (a^c / factorial(c)) * (1/(1-rho)));
Pwait = (a^c / factorial(c)) * (1/(1-rho)) * P0;
Nq_th = (Pwait * rho) / (1 - rho);
Ns_th = Nq_th + a;
Dq_th = Nq_th / lambda;
Ds_th = Dq_th + 1/mu;

% ---------- Plots ----------
figure; plot(MC, Ns_run, 'LineWidth', 1.2); hold on;
yline(Ns_th, 'r--', 'Ns theory');
title('Ns (Average number in system)'); xlabel('MC'); ylabel('Ns');
grid on;

figure; plot(MC, Nq_run, 'LineWidth', 1.2); hold on;
yline(Nq_th, 'r--', 'Nq theory');
title('Nq (Average number in queue)'); xlabel('MC'); ylabel('Nq');
grid on;

figure; plot(completed, Ds_run, 'LineWidth', 1.2); hold on;
yline(Ds_th, 'r--', 'Ds theory');
title('Ds (Average System Delay)'); xlabel('Completed'); ylabel('Ds');
grid on;

figure; plot(completed, Dq_run, 'LineWidth', 1.2); hold on;
yline(Dq_th, 'r--', 'Dq theory');
title('Dq (Average Queue Delay)'); xlabel('Completed'); ylabel('Dq');
grid on;
