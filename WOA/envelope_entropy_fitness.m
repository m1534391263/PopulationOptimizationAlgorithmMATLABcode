function fitness = envelope_entropy_fitness(params, ecg_signal)
    alpha = params(1);
    K = round(params(2)); % ensure K is int
    tau = 0;
    DC = 0;
    init = 1;
    tol = 1e-7;
    
    %% use VMD to denoising
    [u, ~] = VMD(ecg_signal, alpha, tau, K, DC, init, tol);
    denoised_signal = sum(u, 1) - u(1, :);
    
    % calculate the signal envelope
    analytic_signal = hilbert(denoised_signal);
    envelope = abs(analytic_signal);
    
    % permutation entropy of envelope to anew fitness
    fitness = permutation_entropy(envelope, 3, 1);
end
