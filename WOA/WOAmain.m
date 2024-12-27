fprintf('-------------------WOA------------------\n');
% clc;
clear;
close all;

%% input ECG signal
edfFile = ['r01.edf'];
[hdr, data] = edfread(edfFile);
ECG_noisy = data(1, 1:10000);

%% initiate parameters
fs = 1000; % Sampling frequency
t = (0:length(ECG_noisy)-1) / fs; % Time axis

Agents_no = 3;   % number of whales (population size)
Max_iter = 3;    % number of iterations
% dim = 2;          % dimension of the problem (K&a 2paras)
% lb = [4, 500];    % lower bounds (minimum value of K, minimum value of a)
% ub = [10, 2000];  % upper bounds (maximum value of K, maximum value of a)
min_alpha = 500;
max_alpha = 1000;
min_K = 5;
max_K = 10;

% K_range = [4, 10]; % range of K
% alpha_range = [500, 2000]; % range of a
% best_K = 8;
% best_alpha = 796;

tau = 0;
DC = 0;
init = 1;
tol = 1e-7;

%% use WOA to get VMD paras
% [best_alpha, best_K] = woa_optimize_vmd(ECG_noisy, Agents_no, ...
%    min_alpha, max_alpha, min_K, max_K, Max_iter, @fobj);

%% the best paras
% best_K = 8;
% best_alpha = 796;
best_K = 5;
best_alpha = 542.32;
fprintf('best_K = %d, best_alpha = %f\n', best_K, best_alpha);

[u, ~] = VMD(ECG_noisy, best_alpha, tau, best_K, DC, init, tol);
[m,n]=size(u);
% denoised_signal = sum(u,1) - u(1,:);
% denoised_signal = sum(u,1);

%% PCC denoising
% Pearson Correlation Coefficient
% pcc_value = corr(denoising,origin);
% [m,n]=size(u);

denoised_signal = sum(u,1);
for i=2:m
    if corr(u(i,:),ECG_noisy) >= corr(denoised_signal,ECG_noisy)
        denoised_signal = denoised_signal - u(i,:); 
    end
end
% fprintf('---------PCCreconstructed---------\n');
% snr_cc_mse(ECG_noisy,denoised_signal);


%% SNR denoising
% [m,n]=size(u);

denoised_signal = sum(u,1);
for i=1:m
    if snr(u(i,:)) >= snr(denoised_signal)
        denoised_signal = sum(u,1) - u(i,:); 
    end
end
% fprintf('---------SNRreconstructed---------\n');
% snr_cc_mse(ECG_noisy,denoised_signal);

%% figure
denoised_signal = sum(u,1) - u(1,:);
% denoised_signal = denoised_signal - u(1,:) ;
% figure
% for i=1:m
%     subplot(m,1,i);
%     plot(t,u(i,:),'b-','linewidth',1.5)
%     ylabel(['IMF',num2str(i)],'FontName', 'Times New Roman', 'FontSize', 13);
%     axis tight
%     grid on;
% end
% xlabel('Sampling Time/s','FontName', 'Times New Roman', 'FontSize', 18);

% figure;
% subplot(3, 1, 1);
% plot(t, ECG_noisy);
% title('Origin ECG Signal','FontName', 'Times New Roman', 'FontSize', 18);
% ylabel('Vlotage/mV','FontName', 'Times New Roman', 'FontSize', 18);
% grid on;
% 
% subplot(3, 1, 2);
% plot(t, u(1,:));
% title('Noise Signal','FontName', 'Times New Roman', 'FontSize', 18);
% ylabel('Vlotage/mV','FontName', 'Times New Roman', 'FontSize', 18);
% grid on;
% 
% subplot(3, 1, 3);
% %plot(t, denoised_signal-u(1,:));
% plot(t, denoised_signal);
% title('Denoising ECG Signal','FontName', 'Times New Roman', 'FontSize', 18);
% xlabel('Sampling Time/s','FontName', 'Times New Roman', 'FontSize', 18);
% ylabel('Voltage/mV','FontName', 'Times New Roman', 'FontSize', 18);
% grid on;
WOA_denoised = denoised_signal;
noisy_signal = ECG_noisy;
save('../WOA.mat','WOA_denoised','noisy_signal');
%% SNR CC MSE

snr_cc_mse(ECG_noisy,denoised_signal);

