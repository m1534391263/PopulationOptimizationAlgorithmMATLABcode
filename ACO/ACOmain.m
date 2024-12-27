fprintf('-------------------ACO-------------------\n');
%% ECG信号降噪基于蚁群优化算法
clear;
close all;

%% 导入ECG信号
edfFile = 'r01.edf';
[hdr, ACO_ecg_data] = edfread(edfFile);
noisy_ecg = ACO_ecg_data(1, 1:10000);
Fs = 100; % 采样频率

%% 初始化蚁群算法参数
num_ants = 20;        % 蚂蚁数量
max_iter = 100;        % 最大迭代次数
alpha = 1;            % 信息素重要性
beta = 2;             % 启发函数重要性
rho = 0.5;            % 信息素挥发因子
pheromone_init = 1;   % 信息素初始值
lower_bound = [0.5, 5];   % 带通滤波器低频和高频截止频率的下界
upper_bound = [5, 50];    % 带通滤波器低频和高频截止频率的上界

% 初始化信息素矩阵
pheromone = pheromone_init * ones(num_ants, 2);
fitness = zeros(num_ants, 1); % 存储每只蚂蚁的适应度

%% 蚁群算法优化
best_position = [0, 0];
best_fitness = inf;

for iter = 1:max_iter
    % 每只蚂蚁选择滤波器参数
    positions = lower_bound + (upper_bound - lower_bound) .* rand(num_ants, 2);
    for i = 1:num_ants
        low_cutoff = positions(i, 1);
        high_cutoff = positions(i, 2);
        
        % 设计带通滤波器
        [b, a] = butter(2, [low_cutoff, high_cutoff] / (Fs / 2), 'bandpass');
        filtered_ecg = filter(b, a, noisy_ecg);
        
        % 计算适应度：使用负SNR作为目标（最大化SNR）
        fitness(i) = -snr(filtered_ecg, noisy_ecg);
    end
    
    % 更新全局最优解
    [min_fitness, idx] = min(fitness);
    if min_fitness < best_fitness
        best_fitness = min_fitness;
        best_position = positions(idx, :);
    end
    
    % 信息素更新
    for i = 1:num_ants
        pheromone(i, :) = (1 - rho) * pheromone(i, :) + rho / (fitness(i) + 1e-6);
    end
    
    % 打印当前迭代信息
%     disp(['Iteration: ', num2str(iter), ' Best Fitness: ', num2str(-best_fitness)]);
end

%% 使用最佳滤波器参数进行滤波
low_cutoff = best_position(1);
high_cutoff = best_position(2);
[b, a] = butter(2, [low_cutoff, high_cutoff] / (Fs / 2), 'bandpass');
filtered_ecg_final = filter(b, a, noisy_ecg);

%% 结果可视化
figure;

subplot(2,1,1);
plot(noisy_ecg);
title('Noisy ECG Signal');
xlabel('Samples');
ylabel('Amplitude');

subplot(2,1,2); 
plot(filtered_ecg_final);
title(['ACO Filtered ECG Signal (Low: ', num2str(low_cutoff), ' Hz, High: ', num2str(high_cutoff), ' Hz)']);
xlabel('Samples');
ylabel('Amplitude');

ACO_denoised = filtered_ecg_final;
% 显示SNR对比
SNR_before = snr(noisy_ecg);
SNR_after = snr(filtered_ecg_final);
MSE = sqrt(sum((noisy_ecg - filtered_ecg_final).^2)) ./length(noisy_ecg);
R = corrcoef(noisy_ecg, filtered_ecg_final);
r = R(1, 2);
save('../ACO.mat','ACO_denoised');

disp(['降噪前信噪比 (SNR): ', num2str(SNR_before), ' dB']);
disp(['降噪后信噪比 (SNR): ', num2str(SNR_after), ' dB']);
disp(['降噪前后均方误差为：',num2str(MSE)]);
disp(['降噪前后相关系数为：',num2str(r)]);


