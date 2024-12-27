function entropy = permutation_entropy(signal, m, tau)
    % calculate permutation entropy 计算排列熵
    % signal: input signal
    % m: length
    % tau: time delay
    N = length(signal);
    permutations = perms(1:m);
    num_permutations = size(permutations, 1);
    counts = zeros(num_permutations, 1);
    
    for i = 1:(N - (m - 1) * tau)
        pattern = signal(i:tau:(i + (m - 1) * tau));
        [~, order] = sort(pattern);
        for j = 1:num_permutations
            if isequal(order, permutations(j, :))
                counts(j) = counts(j) + 1;
                break;
            end
        end
    end
    
    prob = counts / sum(counts);
    entropy = -sum(prob .* log(prob + eps)); % avoid log(0) no meaning
end


