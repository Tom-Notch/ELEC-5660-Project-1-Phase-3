function Optimal_path = path_from_A_star(map)
    Optimal_path = [];
    size_map = size(map, 1);

    MAX_X = 10;
    MAX_Y = 10;
    MAX_Z = 10;

    % Define the 3D grid map array. Obstacle=-1, Target = 0, Start=1
    MAP = 2 * (ones(MAX_X, MAX_Y, MAX_Z));

    % Initialize MAP with location of the target
    xval = floor(map(size_map, 1));
    yval = floor(map(size_map, 2));
    zval = floor(map(size_map, 3));

    xTarget = xval;
    yTarget = yval;
    zTarget = zval;
    MAP(xval, yval, zval) = 0;

    % Initialize MAP with location of the obstacle
    for i = 2:size_map - 1
        xval = floor(map(i, 1));
        yval = floor(map(i, 2));
        zval = floor(map(i, 3));
        MAP(xval, yval, zval) = -1;
    end

    % Initialize MAP with location of the start point
    xval = floor(map(1, 1));
    yval = floor(map(1, 2));
    zval = floor(map(1, 3));
    xStart = xval;
    yStart = yval;
    zStart = zval;
    MAP(xval, yval, zval) = 1;

    % Main structure in the A* search
    % =====================================================

    % Container storing nodes to be expanded, along with the f score (f=g+h)
    % Each node's (x,y,z) coordinate and its f score is stored in a row For
    % example, queue = [x1, y1, z1, f1; x2, y2, z2, f2; ...; xn, yn, zn, fn]
    queue = [];

    % Arrays for storing the g score of each node, g score of undiscovered
    % nodes is inf
    g = inf(MAX_X, MAX_Y, MAX_Z);

    % Arrays recording whether a node is expanded (popped from the queue) or
    % not expanded: 1, not expanded: 0
    expanded = zeros(MAX_X, MAX_Y, MAX_Z);

    % Arrays recording the parent of each node
    parents = zeros(MAX_X, MAX_Y, MAX_Z, 3);

    % Start your code here
    % ================================================================

    % Define 6 possible movements in the 3D grid (up, down, left, right,
    % forward, backward)
    movements = [1, 0, 0; -1, 0, 0; 0, 1, 0; 0, -1, 0; 0, 0, 1; 0, 0, -1];

    g(xStart, yStart, zStart) = 0;
    queue = [queue; xStart, yStart, zStart, h(xStart, yStart, zStart, xTarget, yTarget, zTarget)];

    while size(queue, 1) > 0
        % Find the node with the smallest f score
        [~, idx] = min(queue(:, 4));
        current_point = [queue(idx, 1), queue(idx, 2), queue(idx, 3)];
        expanded(current_point(1), current_point(2), current_point(3)) = 1;
        queue(idx, :) = [];

        % Check if the target is reached
        if current_point(1) == xTarget && current_point(2) == yTarget && current_point(3) == zTarget
            break;
        end

        % Expand the node
        for move = 1:6

            new_point = current_point + movements(move, :);

            % Check if the new point is within the map
            if new_point(1) < 1 || new_point(1) > MAX_X || new_point(2) < 1 || new_point(2) > MAX_Y || new_point(3) < 1 || new_point(3) > MAX_Z
                continue;
            end

            % Hit an obstacle
            if MAP(new_point(1), new_point(2), new_point(3)) == -1
                continue;
            end

            % Visited
            if expanded(new_point(1), new_point(2), new_point(3)) == 1
                continue;
            end

            % Update the g score and parent only when the resulting path is
            % shorter
            if g(new_point(1), new_point(2), new_point(3)) > g(current_point(1), current_point(2), current_point(3)) + 1
                g(new_point(1), new_point(2), new_point(3)) = g(current_point(1), current_point(2), current_point(3)) + 1;
                queue = [queue; new_point(1), new_point(2), new_point(3), g(new_point(1), new_point(2), new_point(3)) + h(new_point(1), new_point(2), new_point(3), xTarget, yTarget, zTarget)];
                parents(new_point(1), new_point(2), new_point(3), :) = current_point;
            end

        end

    end

    % Reconstruct the path
    current_point = [xTarget, yTarget, zTarget];
    Optimal_path = [current_point];

    while current_point(1) ~= xStart || current_point(2) ~= yStart || current_point(3) ~= zStart
        current_point = squeeze(parents(current_point(1), current_point(2), current_point(3), :))';
        Optimal_path = [current_point; Optimal_path];
    end

    disp("Optimal Path:")
    disp(Optimal_path)

end

function distance = h(x, y, z, xTarget, yTarget, zTarget)
    distance = sqrt((x - xTarget) ^ 2 + (y - yTarget) ^ 2 + (z - zTarget) ^ 2);
end
