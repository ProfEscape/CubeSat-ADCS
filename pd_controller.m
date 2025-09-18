function T_cmd = pd_controller(q_est, omega_body, q_des, params)
% PD_CONTROLLER - Implements quaternion feedback PD control as per ADCS doc
%
% Inputs:
%   q_est     - Estimated quaternion [q0; q1; q2; q3] (scalar-first), body <- orbit
%   omega_body- Angular velocity in body frame [rad/s] (3x1)
%   q_des     - Desired constant quaternion [4x1]
%   params    - Struct with fields:
%                 Kp : proportional gain (scalar or 3x3)
%                 Kd : derivative gain (scalar or 3x3)
%
% Output:
%   T_cmd     - Commanded control torque [N·m] (3x1), to be mapped to magnetorquers

    % --- Error quaternion: δq = q_des^{-1} ⊗ q_est (Eq. 34)
    qd_inv = [q_des(1); -q_des(2:4)];       % inverse of desired quaternion
    delta_q = quatmultiply(qd_inv', q_est'); % Aerospace Toolbox
    delta_q = delta_q(:);

    % Ensure positive scalar part sign convention
    signscalar = sign(delta_q(1));
    if signscalar == 0
        signscalar = 1;
    end

    % Vector part of error quaternion
    delta_qv = delta_q(2:4);

    % --- PD control law (Eq. 35)
    T_cmd = - signscalar * (params.Kp * delta_qv) ...
            - (params.Kd * omega_body);
end
