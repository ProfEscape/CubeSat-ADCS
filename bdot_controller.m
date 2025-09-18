function [m_cmd, tau_m, state] = bdot_controller(B_meas, params, state)
% BDOT_CONTROLLER - Implements B-dot control exactly as per ADCS design doc
%
% Inputs:
%   B_meas  - Current magnetometer measurement [Tesla] (3x1)
%   params  - Struct with fields:
%               Kp_bdot   : controller gain
%               m_max     : max dipole per axis [A·m^2] (3x1)
%   state   - Persistent storage for previous B field
%             state.B_prev : previous magnetometer reading
%             state.dt     : time step [s] (default = 0.1s)
%
% Outputs:
%   m_cmd   - commanded dipole vector [A·m^2]
%   tau_m   - magnetic torque [N·m]
%   state   - updated state (stores B_prev for next call)

    if nargin < 3 || isempty(state)
        % initialize state if first run
        state.B_prev = B_meas;
        state.dt = 0.1; % [s] as given in design doc
    end

    % --- Equation (32): bdot = (b2 - b1)/dt ---
    Bdot = (B_meas - state.B_prev) / state.dt;

    % --- Equation (31): m = -K * bdot ---
    m_cmd = -params.Kp_bdot .* Bdot;

    % --- Algorithm 6: dipole scaling ---
    for i = 1:3
        if m_cmd(i) > params.m_max(i)
            m_cmd(i) = params.m_max(i);
        elseif m_cmd(i) < -params.m_max(i)
            m_cmd(i) = -params.m_max(i);
        end
    end

    % --- Magnetic torque: tau_m = m × b ---
    tau_m = cross(m_cmd, B_meas);

    % --- Update state ---
    state.B_prev = B_meas;
end

