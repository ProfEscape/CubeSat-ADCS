function tau_aero = aerodynamic_torque(rho, v_rel_body, Cd, A_proj, r_cp, r_cg)
% AERODYNAMIC_TORQUE
% Compute aerodynamic disturbance torque
%
% Inputs:
%   rho        - atmospheric density [kg/m^3]
%   v_rel_body - relative velocity in body frame [m/s, 3x1]
%   Cd         - drag coefficient (dimensionless)
%   A_proj     - projected area normal to v_rel [m^2]
%   r_cp       - center of pressure vector in body frame [m]
%   r_cg       - center of gravity vector in body frame [m]
%
% Output:
%   tau_aero   - 3x1 aerodynamic torque [N*m]

v_rel = v_rel_body(:);
v_norm = norm(v_rel);
if v_norm < 1e-6
    tau_aero = zeros(3,1);
    return;
end

F_aero = -0.5 * rho * v_norm^2 * Cd * A_proj * (v_rel/v_norm);
r_arm  = r_cp - r_cg;
tau_aero = cross(r_arm, F_aero);
end
