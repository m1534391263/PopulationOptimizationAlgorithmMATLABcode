tic;
filename = 'Denoise Functions/WOA/WOAmain.m';
run(filename);
WOA_endtime = toc;
disp(['runtime:',num2str(WOA_endtime)]);

tic;
run('Denoise Functions/ABC/ABCmain.m');
ABC_endtime = toc;
disp(['runtime:',num2str(ABC_endtime)]);

tic;
run('Denoise Functions/ARO/AROmain.m');
ARO_endtime = toc;
disp(['runtime:',num2str(ARO_endtime)]);

tic;
run('Denoise Functions/ACO/ACOmain.m');
ACO_endtime = toc;
disp(['runtime:',num2str(ACO_endtime)]);

load('Denoise Functions/WOA.mat');
load('Denoise Functions/ABC.mat');
load('Denoise Functions/ARO.mat');
load('Denoise Functions/ACO.mat');


subplot(5, 1, 1);
plot(noisy_signal);
ylabel('Original Signal','FontName', 'Times New Roman', 'FontSize', 18);

grid on;
subplot(5, 1, 2);
plot(WOA_denoised);
ylabel('WOA','FontName', 'Times New Roman', 'FontSize', 18);
grid on;
subplot(5, 1, 3);
plot(ABC_denoised);
ylabel('Amplitude','FontName', 'Times New Roman', 'FontSize', 28);
grid on;
subplot(5, 1, 4);
plot(ARO_denoised);
ylabel('ARO','FontName', 'Times New Roman', 'FontSize', 18);
grid on;
subplot(5, 1, 5);
plot(ACO_denoised);
ylabel('ACO','FontName', 'Times New Roman', 'FontSize', 18);
grid on;
han = axes('Visible', 'off'); % 创建一个不可见的全局坐标轴
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';
title('Denoising ECG Signal','FontName', 'Times New Roman', 'FontSize', 28);
xlabel('Time (s)','FontName', 'Times New Roman', 'FontSize', 28); % 通用横坐标标签
ylabel('ABC','FontName', 'Times New Roman', 'FontSize', 18); % 通用纵坐标标签