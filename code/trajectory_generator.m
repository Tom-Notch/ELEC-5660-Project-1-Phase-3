function s_des = trajectory_generator(t, path, h, map)
    % https://www.mathworks.com/help/nav/ref/waypointtrajectory-system-object.
    % html

    time_tol = 25;
    poly_degree = 7; %! polynomial degree, 7 for minimum snap

    persistent times t_variables coefficient_variables segment_functions coefficients

    if nargin > 1 % pre-process can be done here (given waypoints). Pre-define the entire trajectory.
        num_segments = size(path, 1) - 1;

        % Initialize variables
        t_variables = cell(3, 1);
        segment_functions = cell(num_segments, 1);
        coefficient_variables = cell(num_segments, 1);
        coefficients = cell(num_segments, 1);

        times = allocate_time(path, time_tol); % time_tol = 25

        for i_dimension = 1:size(path, 2)
            [t_variables{i_dimension}, coefficient_variables{i_dimension}, segment_functions{i_dimension}, coefficients{i_dimension}] = solve_trajectory(path(:, i_dimension), times, poly_degree);
        end

        % visualize the 2D grid map
        subplot(h);
        % start point
        plot3(map(1, 1) - 0.5, map(1, 2) - 0.5, map(1, 3) - 0.5, 'k.');
        hold on;
        % obstacles
        for obs_cnt = 2:size(map, 1) - 1
            plot3([map(obs_cnt, 1) - 0.2 map(obs_cnt, 1) - 0.8], [map(obs_cnt, 2) - 0.2 map(obs_cnt, 2) - 0.8], [map(obs_cnt, 3) map(obs_cnt, 3)], 'k-');
            hold on;
            plot3([map(obs_cnt, 1) - 0.2 map(obs_cnt, 1) - 0.8], [map(obs_cnt, 2) - 0.8 map(obs_cnt, 2) - 0.2], [map(obs_cnt, 3) map(obs_cnt, 3)], 'k-');
            hold on;
            ox = map(obs_cnt, 1) - 0.9;
            oy = map(obs_cnt, 2) - 0.9;
            oz = map(obs_cnt, 3) - 0.9;
            plotcube([0.8, 0.8, 0.8], [ox, oy, oz], 1, [0.7, 0.7, 0.7]);
            grid minor
            set(gca, 'xtick', [-100:1:100])
            set(gca, 'ytick', [-100:1:100])
            grid off;
            grid on;
            axis equal;
            axis ([-1 6 -1 10 0 4]);
            hold on;
        end

        % target point
        plot3(map(obs_cnt + 1, 1) - 0.5, map(obs_cnt + 1, 2) - 0.5, map(obs_cnt + 1, 3) - 0.5, 'r*');
        hold on;

    else % output desired trajectory here (given time)
        s_des = zeros(11, 1);

        % Determine which segment
        i_segment = find(t >= times(1:end - 1) & t < times(2:end), 1, 'last');

        if isempty(i_segment) % If t is exactly the last time, use the last segment
            i_segment = length(times) - 1;
        end

        for i_dimension = 1:3
            position_function = sym(segment_functions{i_dimension}{i_segment});
            velocity_function = sym(diff(segment_functions{i_dimension}{i_segment}, t_variables{i_dimension}, 1));
            acceleration_function = sym(diff(segment_functions{i_dimension}{i_segment}, t_variables{i_dimension}, 2));

            s_des(i_dimension) = double(subs(position_function, [t_variables{i_dimension}, coefficient_variables{i_dimension}(i_segment, :)], [t, coefficients{i_dimension}{i_segment}]));
            s_des(i_dimension + 3) = double(subs(velocity_function, [t_variables{i_dimension}, coefficient_variables{i_dimension}(i_segment, :)], [t, coefficients{i_dimension}{i_segment}]));
            s_des(i_dimension + 6) = double(subs(acceleration_function, [t_variables{i_dimension}, coefficient_variables{i_dimension}(i_segment, :)], [t, coefficients{i_dimension}{i_segment}]));
        end

    end

end

function times = allocate_time(path, total_time)
    % Calculate the Euclidean distances between consecutive waypoints
    segment_distances = sqrt(sum(diff(path) .^ 2, 2));
    total_distance = sum(segment_distances);

    % Calculate the fraction of total time for each segment based on its
    % length
    time_ratios = segment_distances / total_distance;

    % Allocate time to each segment
    segment_times = total_time * time_ratios;

    % Calculate the cumulative time at each waypoint
    times = [0; cumsum(segment_times)];
end

function [t, coefficient_variables, segment_functions, coefficients] = solve_trajectory(path, times, poly_degree)
    num_waypoints = size(path, 1);
    num_segments = num_waypoints - 1;

    long_d_dimension = (poly_degree + 1) * num_segments;

    segment_functions = cell(num_segments, 1);
    Q = cell(num_segments, 1);
    Q_matrix = zeros(long_d_dimension, long_d_dimension);
    A = cell(num_waypoints, 1);
    A_matrix = zeros(long_d_dimension, long_d_dimension);

    syms t real;
    coefficient_variables = sym('p', [num_segments, poly_degree + 1]);

    for i = 1:num_segments
        % Define the symbolic coefficient variables for the i-th segment
        segment_coefficients = coefficient_variables(i, :);
        poly = sum(segment_coefficients .* (t - times(i)) .^ (0:poly_degree));

        segment_functions{i} = matlabFunction(poly, 'Vars', [t, segment_coefficients]);
        cost_function = int(diff(segment_functions{i}, t, (poly_degree + 1) / 2) ^ 2, t, times(i), times(i + 1));

        % Initialize the Q matrix for the segment
        Q{i} = sym(zeros(poly_degree + 1));

        % Extract the terms and coefficients from the integrated cost
        % expression
        [coeff_terms, vars] = coeffs(cost_function, segment_coefficients);

        for j = (poly_degree + 1) / 2 + 1:poly_degree + 1

            for k = j:poly_degree + 1
                index = find(ismember(vars, segment_coefficients(k) * segment_coefficients(j)));

                if j == k
                    Q{i}(j, k) = coeff_terms(index);
                else
                    Q{i}(j, k) = coeff_terms(index) / 2;
                    Q{i}(k, j) = coeff_terms(index) / 2;
                end

            end

        end

        % Fill in Q matrix
        Q_matrix((i - 1) * (poly_degree + 1) + 1:i * (poly_degree + 1), (i - 1) * (poly_degree + 1) + 1:i * (poly_degree + 1)) = Q{i};

        % Construct the A matrix for the segment
        A{i} = zeros(poly_degree + 1, poly_degree + 1);

        for constraint_degree = 0:(poly_degree - 1) / 2
            func = sym(diff(segment_functions{i}, t, constraint_degree));
            % start constraint
            start_constraint = subs(func, t, times(i));
            [start_coefficient, vars] = coeffs(start_constraint, segment_coefficients);

            for j = 1:length(segment_coefficients)
                % Find the index of the current coefficient in the vars array
                var_idx = find(vars == segment_coefficients(j), 1);

                % If the coefficient is found, place it in the A matrix row
                if ~isempty(var_idx)
                    A{i}(constraint_degree + 1, j) = start_coefficient(var_idx);
                end

            end

            % end constraint
            end_constraint = subs(func, t, times(i + 1));
            [end_coefficient, vars] = coeffs(end_constraint, segment_coefficients);

            for j = 1:length(segment_coefficients)
                % Find the index of the current coefficient in the vars array
                var_idx = find(vars == segment_coefficients(j), 1);

                % If the coefficient is found, place it in the A matrix row
                if ~isempty(var_idx)
                    A{i}((poly_degree + 1) / 2 + constraint_degree + 1, j) = end_coefficient(var_idx);
                end

            end

        end

        % Fill in A matrix
        A_matrix((i - 1) * (poly_degree + 1) + 1:i * (poly_degree + 1), (i - 1) * (poly_degree + 1) + 1:i * (poly_degree + 1)) = A{i};

    end

    % Construct C matrix
    C = zeros(long_d_dimension, long_d_dimension - (num_segments - 1) * (poly_degree + 1) / 2);

    % set the very start to 1
    C(1, 1) = 1;

    for i = 1:num_segments - 1
        C((poly_degree + 1) * i - (poly_degree - 1) / 2, i + 1) = 1;
        C((poly_degree + 1) * i + 1, i + 1) = 1;
    end

    % end point 3D location
    C(long_d_dimension - (poly_degree - 1) / 2, num_segments + 1) = 1;

    i_set_index = num_segments + 1;

    % start point derivatives
    for i = 2:(poly_degree + 1) / 2
        C(i, i_set_index + i - 1) = 1;
    end

    i_set_index = i_set_index + (poly_degree - 1) / 2;

    % end point derivatives
    for i = 2:(poly_degree + 1) / 2
        C(long_d_dimension - (poly_degree + 1) / 2 + i, i_set_index + i - 1) = 1;
    end

    i_set_index = i_set_index + (poly_degree - 1) / 2;

    for i_segment = 1:num_segments - 1

        for i_degree = 2:(poly_degree + 1) / 2
            C((i_segment - 1) * (poly_degree + 1) + (poly_degree + 1) / 2 + i_degree, i_set_index + (i_segment - 1) * 3 + i_degree - 1) = 1;
            C(i_segment * (poly_degree + 1) + i_degree, i_set_index + (i_segment - 1) * 3 + i_degree - 1) = 1;
        end

    end

    % Construct d_c vector
    d_c = zeros(i_set_index, 1); % Implicitly setting start and end derivatives to 0
    d_c(1:num_segments + 1) = path;

    % Construct R matrix
    R = C' * inv(A_matrix') * Q_matrix * inv(A_matrix) * C;

    R_uu = R(i_set_index + 1:end, i_set_index + 1:end);
    R_cu = R(1:i_set_index, i_set_index + 1:end);

    % Solve for d_u
    d_u = -inv(R_uu) * R_cu' * d_c;

    d = [d_c; d_u];

    % Translate to coefficients
    p = inv(A_matrix) * C * d;

    coefficients = cell(num_segments, 1);

    for i_num_segments = 1:num_segments
        coefficients{i_num_segments} = p((i_num_segments - 1) * (poly_degree + 1) + 1:i_num_segments * (poly_degree + 1))';
    end

end
