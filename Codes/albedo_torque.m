function tau_alb = albedo_torque(Cr, A_proj, r_cp, r_cg, sun_vec_body, albedo_coeff)
% ALBEDO_TORQUE
% Compute Earth albedo torque (simplified)
%
% Inputs:
%   Cr           - reflectivity coefficient
%   A_proj       - projected area [m^2]
%   r_cp, r_cg   - position vectors [m]
%   sun_vec_body - Sun vector in body frame (approx. direction of incoming flux)
%   albedo_coeff - fraction of solar pressure due to Earth albedo (~0.3 typical)
%
% Output:
%   tau_alb      - 3x1 torque [N*m]

P0 = 4.56e-6; % N/m^2

F_alb = -P0 * albedo_coeff * Cr * A_proj * sun_vec_body(:);
r_arm = r_cp - r_cg;
tau_alb = cross(r_arm, F_alb);
end
