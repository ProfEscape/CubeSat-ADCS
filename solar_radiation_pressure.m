function tau_srp = solar_radiation_pressure(Cr, A_proj, r_cp, r_cg, sun_vec_body, eclipse_flag)
% SOLAR_RADIATION_PRESSURE
% Compute SRP disturbance torque
%
% Inputs:
%   Cr          - reflectivity coefficient (1 ~ 2)
%   A_proj      - projected area facing Sun [m^2]
%   r_cp        - center of pressure [m]
%   r_cg        - center of gravity [m]
%   sun_vec_body- 3x1 unit vector from sat to Sun in body frame
%   eclipse_flag- 1 if in eclipse, 0 if sunlit
%
% Output:
%   tau_srp     - 3x1 SRP torque [N*m]

P0 = 4.56e-6; % N/m^2 @ 1 AU

if eclipse_flag
    tau_srp = zeros(3,1);
    return;
end

F_srp = -P0 * Cr * A_proj * sun_vec_body(:);
r_arm = r_cp - r_cg;
tau_srp = cross(r_arm, F_srp);
end
