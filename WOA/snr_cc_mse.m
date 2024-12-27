function [ ] = snr_cc_mse(ECG_noisy,denoised_signal)
%% SNR MSE CC
    SNR_before = snr(ECG_noisy);
    SNR_after = snr(denoised_signal);
    MSE = sqrt(sum((ECG_noisy - denoised_signal).^2)) ./length(ECG_noisy);
    R = corrcoef(ECG_noisy, denoised_signal);
    r = R(1, 2);
    
    
    disp(['SNR_before: ', num2str(SNR_before), ' dB']);
    disp(['SNR_after: ', num2str(SNR_after), ' dB']);
    disp(['MSE：',num2str(MSE)]);
    disp(['CC：',num2str(r)]);
end

