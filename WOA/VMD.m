%% STANDARD
function [u, u_hat, omega] = VMD(signal, alpha, tau, K, DC, init, tol)
% Input and Parameters:
% ---------------------
% signal  - the time domain signal (1D) to be decomposed
% alpha   - the balancing parameter of the data-fidelity constraint
% tau     - time-step of the dual ascent ( pick 0 for noise-slack )
% K       - the number of modes to be recovered
% DC      - true if the first mode is put and kept at DC (0-freq)
% init    - 0 = all omegas start at 0
%                    1 = all omegas start uniformly distributed
%                    2 = all omegas initialized randomly
% tol     - tolerance of convergence criterion; typically around 1e-6
%
% Output:
% -------
% u       - the collection of decomposed modes
% u_hat   - spectra of the modes
% omega   - estimated mode center-frequencies

% Period and sampling frequency of input signal
save_T = length(signal);
fs = 1/save_T;

% extend the signal by mirroring
T = save_T;
f_mirror(1:T/2) = signal(T/2:-1:1);
f_mirror(T/2+1:3*T/2) = signal;
f_mirror(3*T/2+1:2*T) = signal(T:-1:T/2+1);
f = f_mirror;

% Time Domain 0 to T (of mirrored signal)
T = length(f);
t = (1:T)/T;

% Spectral Domain discretization
freqs = t-0.5-1/T;

% Maximum number of iterations (if not converged yet, then it won't anyway)
N = 500;

% For future generalizations: individual alpha for each mode
Alpha = alpha*ones(1,K);

% Construct and center f_hat
f_hat = fftshift((fft(f)));
f_hat_plus = f_hat;
f_hat_plus(1:T/2) = 0;

% matrix keeping track of every iterant // could be discarded for mem
u_hat_plus = zeros(N, length(freqs), K);

% Initialization of omega_k
omega_plus = zeros(N, K);
switch init
    case 1
        for i = 1:K
            omega_plus(1,i) = (0.5/K)*(i-1);
        end
    case 2
        omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
    otherwise
        omega_plus(1,:) = 0;
end

% if DC mode imposed, set its omega to 0
if DC
    omega_plus(1,1) = 0;
end

% start with empty dual variables
lambda_hat = zeros(N, length(freqs));

% other inits
uDiff = tol+eps; % update step
n = 1; % loop counter
sum_uk = 0; % accumulator



% ----------- Main loop for iterative updates




while ( uDiff > tol &&  n < N ) % not converged and below iterations limit
    
    % update first mode accumulator
    k = 1;
    sum_uk = u_hat_plus(n,:,K) + sum_uk - u_hat_plus(n,:,1);
    
    % update spectrum of first mode through Wiener filter of residuals
    u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
    
    % update first omega if not held at 0
    if ~DC
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
    end
    
    % update of any other mode
    for k=2:K
        
        % accumulator
        sum_uk = u_hat_plus(n+1,:,k-1) + sum_uk - u_hat_plus(n,:,k);
        
        % mode spectrum
        u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
        
        % center frequencies
        mega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
        
    end
    
    % Dual ascent
    lambda_hat(n+1,:) = lambda_hat(n,:) + tau*(sum(u_hat_plus(n+1,:,:),3) - f_hat_plus);
    
    % loop counter
    n = n+1;
    
    % converged yet?
    uDiff = eps;
    for i=1:K
        uDiff = uDiff + 1/T*(u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i))*conj((u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i)))';
    end
    uDiff = abs(uDiff);
    
end


%------ Postprocessing and cleanup


% discard empty space if converged early
N = min(N,n);
omega = omega_plus(1:N,:);

% Signal reconstruction
u_hat = zeros(T, K);
u_hat((T/2+1):T,:) = squeeze(u_hat_plus(N,(T/2+1):T,:));
u_hat((T/2+1):-1:2,:) = squeeze(conj(u_hat_plus(N,(T/2+1):T,:)));
u_hat(1,:) = conj(u_hat(end,:));

u = zeros(K,length(t));

for k = 1:K
    u(k,:)=real(ifft(ifftshift(u_hat(:,k))));
end

% remove mirror part
u = u(:,T/4+1:3*T/4);

% recompute spectrum
clear u_hat;
for k = 1:K
    u_hat(:,k)=fftshift(fft(u(k,:)))';
end

end

%% OTHER
% function [u, omega] = VMD(f, alpha, tau, K, DC, init, tol)
%     % VMD核心函数
% 
%     % 默认参数和信号扩展
%     save_T = length(f);
%     T = save_T * 2; % 确保足够长度避免越界
%     if save_T < 2
%         error('信号长度过短，无法执行VMD分解');
%     end
% 
%     % 扩展信号（避免越界）
%     f_mirror = [f(end:-1:1), f, f(end:-1:1)];
%     f = f_mirror;
% 
%     % 初始化
%     freqs = (1:T) / T - 0.5 - 1 / T;
%     f_hat = fftshift(fft(f));
%     f_hat_plus = f_hat;
%     f_hat_plus(1:T/2) = 0;
% 
%     u_hat_plus = zeros(500, T, K);
%     omega_plus = zeros(500, K);
%     omega_plus(1, :) = linspace(0, 0.5, K); % 初始频率分布
% 
%     lambda_hat = zeros(500, T);
% 
%     % 主循环
%     for n = 1:500
%         for k = 1:K
%             % 更新模态频率
%             sum_uk = sum(u_hat_plus(n, :, :), 3) - u_hat_plus(n, :, k);
%             u_hat_plus(n+1, :, k) = ...
%                 (f_hat_plus - sum_uk - lambda_hat(n, :) / 2) ./ ...
%                 (1 + alpha * (freqs - omega_plus(n, k)).^2);
% 
%             omega_plus(n+1, k) = ...
%                 (freqs * abs(u_hat_plus(n+1, :, k)').^2) / ...
%                 sum(abs(u_hat_plus(n+1, :, k)').^2);
%         end
% 
%         % 双重上升
%         lambda_hat(n+1, :) = lambda_hat(n, :) + tau * ...
%             (sum(u_hat_plus(n+1, :, :), 3) - f_hat_plus);
% 
%         % 收敛检查
%         if max(abs(lambda_hat(n+1, :) - lambda_hat(n, :))) < tol
%             break;
%         end
%     end
% 
%     % 逆傅里叶重构信号
%     u_hat = ifft(ifftshift(u_hat_plus(n, :, :), 2), [], 2);
%     u = real(u_hat(:, save_T/2+1:3*save_T/2));
%     omega = omega_plus(1:n, :);
% end
