function [best_alpha, best_K] = woa_optimize_vmd(ecg_signal, num_agents, ...
    min_alpha, max_alpha, min_K, max_K, Max_iter, fobj)
    % init whales position(para)
    agents = rand(num_agents, 2); % 1st_alpha,2nd_K range(0-1)
    % zoom alpha to[min_alpha, max_alpha]
    agents(:, 1) = agents(:, 1) * (max_alpha - min_alpha) + min_alpha; 
    % zoom K to[min_K, max_K]
    agents(:, 2) = agents(:, 2) * (max_K - min_K) + min_K; 
    
    % init bound
    lb = [min_alpha, min_K]; % paras lower bound
    ub = [max_alpha, max_K]; % paras upper bound
    
    % init the best solution
    best_fitness = inf; % init fitness is infinity
    best_solution = [];
    
    % calculate the best fitness
    for i = 1:num_agents
        fitness = fobj(agents(i, :), ecg_signal); 
        % use fobj to calculate fitness
        if fitness < best_fitness
            best_fitness = fitness;
            best_solution = agents(i, :);
        end
    end
    
    % WOA main loop
    for t = 1:Max_iter
        % anew a&a2
        a = 2 - t * (2 / Max_iter);
        a2 = -1 + t * (-1 / Max_iter);
        
        for i = 1:num_agents
            % anew position
            r1 = rand();
            r2 = rand();
            A = 2 * a * r1 - a;
            C = 2 * r2;
            p = rand();
            D_Leader = abs(C * best_solution - agents(i, :));
            
            if p < 0.5
                if abs(A) >= 1
                    % random choose a whale to anew position
                    rand_leader_index = randi([1, num_agents]);
                    X_rand = agents(rand_leader_index, :);
                    D_X_rand = abs(C * X_rand - agents(i, :));
                    new_position = X_rand - A * D_X_rand;
                else
                    % Referring to the current global optimal solution
                    new_position = best_solution - A * D_Leader;
                end
            else
                % Adopting a spiral update approach
                l = (a2 - 1) * rand + 1;
                b = 1;
                new_position = D_Leader .* exp(b .* l) .* cos(l .* 2 * pi) + best_solution;
            end
            
            % to ensure new_position in boundary
            new_position(1) = max(min(new_position(1), max_alpha), min_alpha); % alpha边界处理
            new_position(2) = round(max(min(new_position(2), max_K), min_K)); % K边界处理，且确保为整数
            
            % calculate new_fitness
            new_fitness = fobj(new_position, ecg_signal);
            
            % anew whale position
            current_fitness = fobj(agents(i, :), ecg_signal);
            if new_fitness < current_fitness
                agents(i, :) = new_position;
            end
            
            % anew best_solution
            if new_fitness < best_fitness
                best_fitness = new_fitness;
                best_solution = new_position;
            end
        end
        
        % output iteration infos
        fprintf(['Iteration_number_t = %d, best_alpha = %.4f, best_K = %d,' ...
            ' optimal_fitness = %.6f\n'], ...
                t, best_solution(1), round(best_solution(2)), best_fitness);
    end
    
    % return best_alpha & best_K
    best_alpha = best_solution(1);
    best_K = round(best_solution(2));
end
