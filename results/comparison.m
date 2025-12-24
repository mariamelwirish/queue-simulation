clear; clc; close all;

% ALGORITHMS
files = [
    struct('name','fifo','file','mmc_fifo.csv')
    struct('name','edf','file','mmc_edf.csv')
    struct('name','sjf','file','mmc_sjf.csv')
];

basePath = "csv/";  

%  LOAD ALL DATA
numAlgs = length(files);
data = cell(numAlgs,1);

for i = 1:numAlgs
    fullPath = fullfile(basePath, files(i).file);
    data{i} = readtable(fullPath);

    % Basic sanity check
    requiredCols = {'rho','Ns_sim','Nq_sim','Ds_sim','Dq_sim','pExpired_sim'};
    if ~all(ismember(requiredCols, data{i}.Properties.VariableNames))
        error("CSV %s does not have the required structure.", files(i).file);
    end
end

rho = data{1}.rho;   % assume same rho for all

%  METRICS TO COMPARE
metrics = {
    'Ns_sim',        'Average number in system (N_s)';
    'Nq_sim',        'Average number in queue (N_q)';
    'Ds_sim',        'Average system delay (D_s)';
    'Dq_sim',        'Average queue delay (D_q)';
    'pExpired_sim',  'Probability of expiration';
    'pServed', 'Probability of being served';
};

%  PLOTTING
for m = 1:size(metrics,1)

    metricCol  = metrics{m,1};
    metricName = metrics{m,2};

    figure;
    hold on; grid on;

    for i = 1:numAlgs
        plot(rho, data{i}.(metricCol), ...
            '-o', ...                 % line + circle markers
            'LineWidth', 2, ...
            'MarkerSize', 6, ...
            'DisplayName', files(i).name);

    end

    xlabel('\rho');
    ylabel(metricName);
    title([metricName]);
    legend('Location','best');
end
