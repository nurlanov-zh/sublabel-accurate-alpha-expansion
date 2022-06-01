function [labels_matrix, quadratic, count_non_convex ...
    ] = local_quadratic_data_approx(label_inds_gco, label_space_double, ...
    data_cost, label_neighbor_size, bound_for_adaptive)
% Approximate data term locally by convex quadratic functions

    [L, N] = size(data_cost);
    
    left_inds = double(min(max(1, label_inds_gco - 1), ...
                            L - label_neighbor_size + 1));
    right_inds = double(left_inds + label_neighbor_size - 1);
    center_inds = left_inds + 1;
    
    %%%%%%%%%%%%%%%%% Local Quadratic interpolation of data term in 1D case %%%%%%
    % quadratic - for each vertex save coefficients (a,b,c): x*x*a + x*b + c
    
    quadratic = zeros(N, 3);
    count_non_convex = 0;
    for cite=1:N
        move_right = true;
        prev_left = 0;
        prev_right = 0;
        while true
            % create X, Y for approximating data_cost
            num_points_for_approx = right_inds(cite) - left_inds(cite) + 1;
            
            y_costs = zeros(num_points_for_approx, 1);
            x_inds = zeros(num_points_for_approx, 1);
            x_labels = zeros(length(x_inds), 1);
            for i = 1:num_points_for_approx
                x_inds(i) = left_inds(cite) + i - 1;
                y_costs(i) = data_cost(x_inds(i), cite);
                x_labels(i) = label_space_double(x_inds(i));
            end
            
            % normalize X
            % x_labels = (x_labels - 1) / (L - 1);
            % Fit quadratic
            [quadr_coefs, residual] = fit_quadratic(x_labels, y_costs);
            % Get rid of -0 as a coefficient
            if abs(quadr_coefs(1)) <= 1e-8
                quadr_coefs(1) = 0;
            end
            
            % If convex, continue loop
            if quadr_coefs(1) >= 0 && residual < 0.05
                % Save previous iteration
                quadratic(cite, :) = quadr_coefs;
                if left_inds(cite) == 1 && right_inds(cite) == L
                    break;
                end
                % Bound the neighborhood
                if num_points_for_approx >= bound_for_adaptive
                    break
                end
                prev_left = left_inds(cite);
                prev_right = right_inds(cite);
                % Move borders
                if move_right && right_inds(cite) < L
                    right_inds(cite) = min(L, right_inds(cite) + 1);
                    if left_inds(cite) > 1
                        move_right = ~move_right;
                    end
                else
                    left_inds(cite) = max(1, left_inds(cite) - 1);
                    if right_inds(cite) < L
                        move_right = ~move_right;
                    end
                end
            else
                % If non-convex
                if prev_left ~= 0 && prev_right ~= 0
                    left_inds(cite) = prev_left;
                    right_inds(cite) = prev_right;
                    break
                end
                count_non_convex = count_non_convex + 1;
                %%% least squares
                % lin_coefs = fit_linear(x_labels, y_costs);
                
                %%% contains original point (mid) and minimal 
                left = 1;
                mid = 2;
                right = 3;
                
                if left_inds(cite) == 1
                    lin_coefs = fit_linear([x_labels(left); x_labels(mid)], ...
                        [y_costs(left); y_costs(mid)]);
                    % right_inds(cite) = center_inds(cite);
                elseif right_inds(cite) == L
                    lin_coefs = fit_linear([x_labels(mid); x_labels(right)], ...
                        [y_costs(mid); y_costs(right)]);
                    % left_inds(cite) = center_inds(cite);
                else
                    if y_costs(left) <= y_costs(right)
                        lin_coefs = fit_linear([x_labels(left); x_labels(mid)], ...
                            [y_costs(left); y_costs(mid)]);
                        % right_inds(cite) = center_inds(cite);
                    else
                        lin_coefs = fit_linear([x_labels(mid); x_labels(right)], ...
                            [y_costs(mid); y_costs(right)]);
                        % left_inds(cite) = center_inds(cite);
                    end
                end
                
                quadratic(cite, :) = [0; lin_coefs];
                break
            end
        end
    end
    
    labels_matrix = [left_inds center_inds right_inds];
end

