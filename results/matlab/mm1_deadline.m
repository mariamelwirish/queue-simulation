lambda = 0.18;     % arrival rate
mu     = 0.25;     % service rate
delta  = 0.12;     % deadline timer rate (for reference; theory below ignores deadlines)

T = readtable('mm1_timeseries.csv');   % per-event running avgs

MC        = T.MC;
completed = T.completed;
Ns_run    = T.Ns_running;
Nq_run    = T.Nq_running;
Ds_run    = T.Ds_running;
Dq_run    = T.Dq_running;
pe_run    = T.pe_running;     

rho = lambda / mu;
if rho >= 1
    error('Unstable M/M/1: rho >= 1. Choose lambda < mu.');
end

Ns_th = rho / (1 - rho);
Nq_th = (rho^2) / (1 - rho);
Ds_th = 1 / (mu - lambda);
Dq_th = Nq_th / lambda;

figure; plot(MC, Ns_run, 'LineWidth', 1.2); hold on;
yline(Ns_th, 'r--', 'Ns theory');
title('Ns (Average number in system)'); xlabel('MC'); ylabel('Ns'); grid on;

figure; plot(MC, Nq_run, 'LineWidth', 1.2); hold on;
yline(Nq_th, 'r--', 'Nq theory');
title('Nq (Average number in queue)'); xlabel('MC'); ylabel('Nq'); grid on;

figure; plot(completed, Ds_run, 'LineWidth', 1.2); hold on;
yline(Ds_th, 'r--', 'Ds theory');
title('Ds (Average System Delay)'); xlabel('Completed'); ylabel('Ds'); grid on;

figure; plot(completed, Dq_run, 'LineWidth', 1.2); hold on;
yline(Dq_th, 'r--', 'Dq theory');
title('Dq (Average Queue Delay)'); xlabel('Completed'); ylabel('Dq'); grid on;

figure; plot(completed, pe_run, 'LineWidth', 1.2);
title('P(served after expiry) — running estimate'); xlabel('Completed'); ylabel('P_{expired-served}'); grid on;


