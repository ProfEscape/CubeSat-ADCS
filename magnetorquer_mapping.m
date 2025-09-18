function [m_cmd_sat, T_magnetic_effective] = magnetorquer_mapping(T_cmd, B_body, params)
% MAGNETORQUER_MAPPING - Torque allocation per ADCS design document
%
% Inputs:
%   T_cmd   - Commanded control torque [N·m] (3x1)
%   B_body  - Magnetic field in body frame [Tesla] (3x1)
%   params  - Struct with fields:
%               m_max      : [3x1] max dipole per axis [A·m^2]
%               residual_m : [3x1] residual dipole [A·m^2]
%
% Outputs:
%   m_cmd_sat            - Saturated dipole command [A·m^2]
%   T_magnetic_effective - Actual torque produced [N·m]

    if norm(B_body) < 1e-12
        m_cmd_sat = [0;0;0];
        T_magnetic_effective = [0;0;0];
        return;
    end

    % --- Normalize magnetic field ---
    b_hat = B_body / norm(B_body);

    % --- Projection: extract component of T_cmd perpendicular to b_hat ---
    T_magnetic = (skew(b_hat)' * skew(b_hat)) * T_cmd;

    % --- Required dipole to produce T_magnetic ---
    m = -cross(T_magnetic, B_body) / (norm(B_body)^2);

    % --- Subtract known residual dipole ---
    m = m - params.residual_m(:);

    % --- Magnetic dipole scaling (Algorithm 6) ---
    for i = 1:3
        if m(i) > params.m_max(i)
            m(i) = params.m_max(i);
        elseif m(i) < -params.m_max(i)
            m(i) = -params.m_max(i);
        end
    end
    m_cmd_sat = m;

    % --- Effective torque with real magnetic field ---
    T_magnetic_effective = cross(m_cmd_sat + params.residual_m(:), B_body);
end

% --- Skew-symmetric helper ---
function S = skew(v)
S = [  0   -v(3)  v(2);
      v(3)   0   -v(1);
     -v(2) v(1)    0 ];
end

