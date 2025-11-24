clear; clc; close all;

set(0, 'DefaultFigureWindowStyle', 'docked');

csv_fifo = '../csv/mmc_fifo.csv';
csv_edf  = '../csv/mmc_edf.csv';
csv_sjf  = '../csv/mmc_sjf.csv';

Tf = readtable(csv_fifo);
Te = readtable(csv_edf);
Ts = readtable(csv_sjf);


rho = Tf.rho;


Ns_fifo   = Tf.Ns_sim;
Nq_fifo   = Tf.Nq_sim;
Ds_fifo   = Tf.Ds_sim;
Dq_fifo   = Tf.Dq_sim;
pExp_fifo = Tf.pExpired_sim;


Ns_edf   = Te.Ns_sim;
Nq_edf   = Te.Nq_sim;
Ds_edf   = Te.Ds_sim;
Dq_edf   = Te.Dq_sim;
pExp_edf = Te.pExpired_sim;


Ns_sjf   = Ts.Ns_sim;
Nq_sjf   = Ts.Nq_sim;
Ds_sjf   = Ts.Ds_sim;
Dq_sjf   = Ts.Dq_sim;
pExp_sjf = Ts.pExpired_sim;


Ns_th = Tf.Ns_th;
Nq_th = Tf.Nq_th;
Ds_th = Tf.Ds_th;
Dq_th = Tf.Dq_th;


outdir = 'plots';
if ~exist(outdir,'dir')
    mkdir(outdir);
end


make_compare_plot3(rho, Ns_fifo, Ns_edf, Ns_sjf, Ns_th, ...
    'Average Number in System (N_s)', 'N_s', fullfile(outdir, 'Ns'));

make_compare_plot3(rho, Nq_fifo, Nq_edf, Nq_sjf, Nq_th, ...
    'Average Number in Queue (N_q)', 'N_q', fullfile(outdir, 'Nq'));

make_compare_plot3(rho, Ds_fifo, Ds_edf, Ds_sjf, Ds_th, ...
    'Average System Delay (D_s)', 'Time', fullfile(outdir, 'Ds'));

make_compare_plot3(rho, Dq_fifo, Dq_edf, Dq_sjf, Dq_th, ...
    'Average Queueing Delay (D_q)', 'Time', fullfile(outdir, 'Dq'));

make_compare_plot3_no_theory(rho, pExp_fifo, pExp_edf, pExp_sjf, ...
    'P(Expired Task Served)', 'Probability', fullfile(outdir, 'pExpired'));

disp('All plots generated, saved, and left open (docked).');

function make_compare_plot3(x, y_fifo, y_edf, y_sjf, y_th, ttl, ylab, filename)
    f = figure; hold on; grid on;

    plot(x, y_fifo, 'ro-', 'LineWidth', 1.6, 'MarkerSize', 6);
    plot(x, y_edf,  'b*-', 'LineWidth', 1.6, 'MarkerSize', 6);
    plot(x, y_sjf,  'gs-', 'LineWidth', 1.6, 'MarkerSize', 6);
    plot(x, y_th,   'k:',  'LineWidth', 2.0);  % Red dotted line

    legend('FIFO (simulation)','EDF (simulation)','SJF (simulation)', ...
           'Erlang-C (M/M/c theory)','Location','best');

    xlabel('\rho');
    ylabel(ylab);
    title(ttl);

    % Save as PNG and PDF
    saveas(f, filename + ".png");
    saveas(f, filename + ".pdf");

end


function make_compare_plot3_no_theory(x, y_fifo, y_edf, y_sjf, ttl, ylab, filename)
    f = figure; hold on; grid on;

    plot(x, y_fifo, 'ro-', 'LineWidth', 1.6, 'MarkerSize', 6);
    plot(x, y_edf,  'b*-', 'LineWidth', 1.6, 'MarkerSize', 6);
    plot(x, y_sjf,  'gs-', 'LineWidth', 1.6, 'MarkerSize', 6);

    legend('FIFO (simulation)','EDF (simulation)','SJF (simulation)', ...
           'Location','best');

    xlabel('\rho');
    ylabel(ylab);
    title(ttl);

    % Save as PNG and PDF
    saveas(f, filename + ".png");
    saveas(f, filename + ".pdf");

    % Do NOT close(f); so you can still see the plot
end
