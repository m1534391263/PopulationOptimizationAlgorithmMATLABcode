function fitness = fobj(params, ecg_signal)
    alpha = params(1);
    K = round(params(2)); % 确保分解层数为整数
    [u, ~] = VMD(ecg_signal, alpha, 0, K, 0, 1, 1e-7); % 调用VMD函数
    denoised_signal = sum(u, 1) - u(1, :); % 简单去噪：保留所有模态
    fitness = envelope_entropy_fitness(params, denoised_signal); % 使用包络熵计算适应度
end