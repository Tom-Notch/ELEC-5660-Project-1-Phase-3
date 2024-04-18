function [F, M] = controller(t, s, s_des)

    global params

    m = params.mass;
    g = params.grav;
    I = params.I;

    % PID parameters
    position_Kp = 5;
    position_Kd = 10;

    attitude_Kp = 0.001;
    attitude_Kd = 100;

    % current state
    current_position = s(1:3); % current position
    current_x = current_position(1);
    current_y = current_position(2);
    current_z = current_position(3);

    current_velocity = s(4:6); % current velocity
    current_x_velocity = current_velocity(1);
    current_y_velocity = current_velocity(2);
    current_z_velocity = current_velocity(3);

    current_attitude_quaternion = s(7:10); % current attitude quaternion
    current_Rot = QuatToRot(current_attitude_quaternion'); % current rotation matrix
    [current_row, current_pitch, current_yaw] = RotToRPY_ZXY(current_Rot); % current roll, pitch, yaw

    % gravity vector in body frame
    g_body = current_Rot * [0; 0; -g];

    current_angular_velocity = s(11:13); % current body angular velocity
    current_row_rate = current_angular_velocity(1);
    current_pitch_rate = current_angular_velocity(2);
    current_yaw_rate = current_angular_velocity(3);

    % desired state
    desired_position = s_des(1:3); % desired position
    desired_x = desired_position(1);
    desired_y = desired_position(2);
    desired_z = desired_position(3);

    desired_velocity = s_des(4:6); % desired velocity
    desired_x_velocity = desired_velocity(1);
    desired_y_velocity = desired_velocity(2);
    desired_z_velocity = desired_velocity(3);

    desired_acceleration = s_des(7:9); % desired acceleration
    desired_x_acceleration = desired_acceleration(1);
    desired_y_acceleration = desired_acceleration(2);
    desired_z_acceleration = desired_acceleration(3);

    desired_yaw = s_des(10); % desired yaw
    desired_yaw_rate = s_des(11); % desired yaw rate

    % Position control
    x_control = position_PID(desired_x_acceleration, desired_x_velocity, desired_x, current_x_velocity, current_x, position_Kp, position_Kd);
    y_control = position_PID(desired_y_acceleration, desired_y_velocity, desired_y, current_y_velocity, current_y, position_Kp, position_Kd);
    z_control = position_PID(desired_z_acceleration, desired_z_velocity, desired_z, current_z_velocity, current_z, position_Kp, position_Kd);

    F = m * (g + z_control);

    desired_row = 1 / g * (x_control * sin(current_yaw) - y_control * cos(current_yaw));
    desired_pitch = 1 / g * (x_control * cos(current_yaw) + y_control * sin(current_yaw));

    desired_row_rate = normalize_angle(desired_row - current_row) / 0.01; % 0.01 is the control interval
    desired_pitch_rate = normalize_angle(desired_pitch - current_pitch) / 0.01; % 0.01 is the control interval
    desired_yaw_rate = normalize_angle(desired_yaw - current_yaw) / 0.01 + desired_yaw_rate; % 0.01 is the control interval, adding desired_yaw_rate for feedforward

    % Attitude control
    row_control = attitude_PID(desired_row_rate, desired_row, current_row_rate, current_row, attitude_Kp, attitude_Kd);
    pitch_control = attitude_PID(desired_pitch_rate, desired_pitch, current_pitch_rate, current_pitch, attitude_Kp, attitude_Kd);
    yaw_control = attitude_PID(desired_yaw_rate, desired_yaw, current_yaw_rate, current_yaw, attitude_Kp, attitude_Kd);

    M = I * [row_control; pitch_control; yaw_control] + cross(current_angular_velocity, I * current_angular_velocity);

end

function linear_control = position_PID(desired_acceleration, desired_velocity, desired_position, current_velocity, current_position, Kp, Kd)
    linear_control = desired_acceleration + Kd * (desired_velocity - current_velocity) + Kp * (desired_position - current_position);
end

function angular_control = attitude_PID(desired_velocity, desired_Rot, current_velocity, current_Rot, Kp, Kd)
    angular_control = Kp * normalize_angle(desired_Rot - current_Rot) + Kd * (desired_velocity - current_velocity);
end

function normalized_angle = normalize_angle(angle)
    % Normalize the angle to be within the range [0, 2*pi)
    normalized_angle = mod(angle, 2 * pi);

    % Adjust angles to be in the range [-pi, pi)
    if normalized_angle > pi
        normalized_angle = normalized_angle - 2 * pi;
    end

end
