fprintf('-------------------ARO-------------------\n');
% ECG信号降噪基于ARO算法
clear;
close all;

% 导入ECG信号
edfFile = 'r01.edf';
[hdr, data1] = edfread(edfFile);
noisy_ecg = data1(1, 1:10000);
Fs = 100; % 采样频率

% ARO算法参数初始化
num_rabbits = 50;  % 兔子数量
max_iter = 100;    % 最大迭代次数
lower_bound = [0.5, 5];   % 带通滤波器低频和高频截止的下界
upper_bound = [5, 50];    % 带通滤波器低频和高频截止的上界

% 初始化兔子的位置（带通滤波器的低频和高频参数）
rabbits = lower_bound + (upper_bound - lower_bound) .* rand(num_rabbits, 2);
fitness = zeros(num_rabbits, 1);  % 用于存储每个兔子的适应度

% 显示初始信息
disp(['总兔子数量: ', num2str(num_rabbits), ', 最大迭代次数: ', num2str(max_iter)]);
disp(['滤波器频率范围: [', num2str(lower_bound(1)), ' Hz, ', num2str(upper_bound(2)), ' Hz]']);


% 适应度函数：计算经过滤波后的ECG信号与原信号之间的相似度
for iter = 1:max_iter
    for i = 1:num_rabbits
        % 使用当前兔子的位置（滤波器的低频和高频截止频率）进行滤波
        low_cutoff = rabbits(i, 1);
        high_cutoff = rabbits(i, 2);
        
        % 设计带通滤波器
        [b, a] = butter(2, [low_cutoff high_cutoff] / (Fs / 2), 'bandpass');
        filtered_ecg = filter(b, a, noisy_ecg);
        
        % 计算适应度：此处使用信号的相似度作为适应度指标（例如SNR、相关性）
        fitness(i) = -snr(filtered_ecg, noisy_ecg);  % SNR越高，适应度越好（目标是最大化SNR）
    end
    
    % ARO算法核心更新过程
    [best_fitness, best_idx] = min(fitness);  % 选择适应度最好的兔子
    best_position = rabbits(best_idx, :);  % 最好的兔子的位置
    best_snr = -best_fitness; % 当前最佳信噪比
    
    for i = 1:num_rabbits
        % 更新兔子的位置（速度和位置更新）
        rabbits(i, :) = rabbits(i, :) + rand(1, 2) .* (best_position - rabbits(i, :));
        
        % 限制兔子的位置在边界内
        rabbits(i, :) = max(rabbits(i, :), lower_bound);
        rabbits(i, :) = min(rabbits(i, :), upper_bound);
    end
    
    % 输出当前迭代的信息
%     disp(['迭代次数: ', num2str(iter), ...
%           ', 最佳适应度: ', num2str(best_fitness), ...
%           ', 最佳SNR: ', num2str(best_snr), ...
%           ' dB (低频: ', num2str(best_position(1)), ...
%           ' Hz, 高频: ', num2str(best_position(2)), ' Hz)']);
end

% 使用最佳参数（滤波器参数）进行最终的ECG信号滤波
low_cutoff = best_position(1);
high_cutoff = best_position(2);
[b, a] = butter(2, [low_cutoff high_cutoff] / (Fs / 2), 'bandpass');
filtered_ecg_final = filter(b, a, noisy_ecg);

% 绘图展示结果
% figure;
% 
% subplot(2,1,1);
% plot(noisy_ecg);
% title('Noisy ECG Signal');
% xlabel('Samples');
% ylabel('Amplitude');
% 
% subplot(2,1,2);
% plot(filtered_ecg_final);
% title(['ARO Filtered ECG Signal (Low: ', num2str(low_cutoff), ' Hz, High: ', num2str(high_cutoff), ' Hz)']);
% xlabel('Samples');
% ylabel('Amplitude');

% 计算并显示SNR对比
SNR_noise = snr(noisy_ecg);
SNR_filtered = snr(filtered_ecg_final);
ARO_denoised =  filtered_ecg_final;
save('../ARO.mat','ARO_denoised');

disp(['降噪前信噪比 (SNR): ', num2str(SNR_noise), ' dB']);
disp(['降噪后信噪比 (SNR): ', num2str(SNR_filtered), ' dB']);
disp(['最佳滤波器频率范围: ', num2str(low_cutoff), ' Hz - ', num2str(high_cutoff), ' Hz']);
