function tau_gg = gravity_gradient_torque(I_body, r_eci_km, q)
% GRAVITY_GRADIENT_TORQUE
% Compute gravity gradient torque in body frame
%
% Inputs:
%   I_body      - 3x3 inertia matrix [kg*m^2]
%   r_eci_km    - 3x1 position in ECI frame [km]
%   q           - 1x4 quaternion [qx qy qz qw] (ECI->body, scalar-last)
%
% Output:
%   tau_gg      - 3x1 gravity gradient torque [N*m]

mu = 3.986004418e5; % km^3/s^2

% Convert r_eci to body frame
A_q = quat_to_rotm(q);      % inertial->body rotation
r_b = A_q * r_eci_km(:);    % km

% Convert km to m
r_b_m = r_b * 1e3;

r_norm = norm(r_b_m);
tau_gg = 3*mu*1e9/(r_norm^5) * cross(r_b_m, I_body*r_b_m); % [N*m]
end
