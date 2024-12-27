fprintf('-------------------ABC-------------------\n');
clear;
close all;
%% 导入ECG信号
edfFile = 'r01.edf';
[hdr, ACO_ecg_data] = edfread(edfFile);
noisy_ecg = ACO_ecg_data(1, 1:10000);
fields_hdr = fieldnames(hdr);
Fs = 100;
%% 初始化人工蜂群算法参数
num_bees = 20;        % 蜜蜂数量（包括雇佣蜂和观察蜂）这里假设它们数量相同
max_iter = 100;        % 最大迭代次数
limit = 100;          % 蜜蜂放弃当前食物源的限制次数（如果适应度没有改进）
lower_bound = [0.5, 5];   % 带通滤波器低频和高频截止频率的下界
upper_bound = [5, 45];    % 带通滤波器低频和高频截止频率的上界

% 初始化蜜蜂位置（即滤波器参数）
positions = lower_bound + (upper_bound - lower_bound) .* rand(num_bees, 2);
fitness = zeros(num_bees, 1); % 存储每只蜜蜂的适应度
trial = zeros(num_bees, 1);   % 记录每只蜜蜂在当前食物源的尝试次数

%% 人工蜂群算法优化
best_position = [0, 0];
best_fitness = -inf;
for i = 1:num_bees
        low_cutoff = positions(i, 1);
        high_cutoff = positions(i, 2);
        
for iter = 1:max_iter
    % 评估每只蜜蜂的适应度
    for i = 1:num_bees
        low_cutoff = positions(i, 1);
        high_cutoff = positions(i, 2);
        % 设计带通滤波器
        [b, a] = butter(2, [low_cutoff, high_cutoff] / (Fs / 2),'bandpass');
        filtered_ecg = filter(b, a, noisy_ecg);
        
    end
end
    % 计算适应度：使用SNR作为目标（最大化SNR）
fitness(i) = snr(filtered_ecg, noisy_ecg); % 注意这里去掉了负号，因为我们要最大化SNR
end
    % 更新全局最优解
    [max_fitness, idx] = max(fitness); % 注意这里使用了max来找到最大的SNR
    if max_fitness > best_fitness
        best_fitness = max_fitness;
        best_position = positions(idx, :);
    end
    % 雇佣蜂和观察蜂阶段（这里合并为一个阶段，使用相同的搜索策略）
    for i = 1:num_bees
        % 选择一个不同于i的随机蜜蜂k
        k = randi(num_bees);
        while k == i
            k = randi(num_bees);
        end
        % 生成新位置（基于当前位置和随机选择的其他蜜蜂位置）
        phi = rand;
        new_position = positions(i, :) + phi * (positions(i, :) - positions(k, :));
        new_position = max(min(new_position, upper_bound), lower_bound); % 确保新位置在边界内
        % 在滤波器设计前检查范围
        low_cutoff = positions(i, 1);
        high_cutoff = positions(i, 2);
%         disp(['low_cutoff: ', num2str(low_cutoff), ', high_cutoff: ', num2str(high_cutoff)]);
        % 评估新位置的适应度
        low_cutoff = new_position(1);
        high_cutoff = new_position(2);
        [b, a] = butter(2, [low_cutoff, high_cutoff] / (Fs / 2), 'bandpass');
        new_filtered_ecg = filter(b, a, noisy_ecg);
        new_fitness = snr(new_filtered_ecg, noisy_ecg);
        % 如果新位置更好，则更新蜜蜂位置
        if new_fitness > fitness(i)
            positions(i, :) = new_position;
            fitness(i) = new_fitness;
            trial(i) = 0; % 重置尝试次数
        else
            trial(i) = trial(i) + 1; % 增加尝试次数
        end
    end
    
    % 打印当前迭代信息
    disp(['Iteration: ', num2str(iter), ' Best Fitness (SNR): ', num2str(best_fitness)]);
% 在滤波器设计前检查范围
low_cutoff = positions(i, 1);
high_cutoff = positions(i, 2);
disp(['low_cutoff: ', num2str(low_cutoff), ', high_cutoff: ', num2str(high_cutoff)]);
%% 使用最佳滤波器参数进行滤波
low_cutoff = best_position(1);
high_cutoff = best_position(2);
[b, a] = butter(2, [low_cutoff, high_cutoff] / (Fs / 2), 'bandpass');
filtered_ecg_final = filter(b, a, noisy_ecg);
% 显示SNR对比
SNR_before = snr(noisy_ecg);
SNR_after = snr(filtered_ecg_final);
MSE = sqrt(sum((noisy_ecg - filtered_ecg_final).^2)) ./length(noisy_ecg);
R = corrcoef(noisy_ecg, filtered_ecg_final);
r = R(1, 2);


disp(['降噪前信噪比 (SNR): ', num2str(SNR_before), ' dB']);
disp(['降噪后信噪比 (SNR): ', num2str(SNR_after), ' dB']);
disp(['降噪前后均方误差为：',num2str(MSE)]);
disp(['降噪前后相关系数为：',num2str(r)]);

%% 结果可视化

% figure;
% subplot(2,1,1);
% plot(noisy_ecg);
% title('Noisy ECG Signal');
% xlabel('Samples');
% ylabel('Amplitude');
% subplot(2,1,2);
% plot(filtered_ecg_final);
% title('Filtered ECG Signal');
% xlabel('Samples');
% ylabel('Amplitude');
ABC_denoised = filtered_ecg_final;
save('../ABC.mat','ABC_denoised');