function [omega_dot, omega_next, q_dot] = attitude_dynamics(I, omega, tau_total, dt, integrate_q, q)
% ATTITUDE_DYNAMICS  Rigid-body rotational dynamics
%
% Inputs:
%   I         - 3x3 inertia matrix (body frame) [kg*m^2]
%   omega     - 3x1 angular velocity (body frame) [rad/s]
%   tau_total - 3x1 total applied torque (body frame) [N*m] (includes control + disturbances)
%   dt        - (optional) time-step for integration (s). If empty, no integration.
%   integrate_q - (optional, default false) boolean: compute q_dot if q provided
%   q         - (optional) 1x4 quaternion [qx qy qz qw], needed only when integrate_q true
%
% Outputs:
%   omega_dot - 3x1 derivative of omega (rad/s^2)
%   omega_next- 3x1 next omega after dt (rad/s) if dt provided, else []
%   q_dot     - 1x4 quaternion derivative (if integrate_q true), else []
%
% Dynamics:
%   omega_dot = I^{-1} * (tau_total - omega x (I*omega))  (DDJF eqn style)
%
if nargin < 5
    integrate_q = false;
end
if nargin < 6
    q = [];
end

% ensure column vectors
omega = omega(:);
tau_total = tau_total(:);

% coriolis term
c = cross(omega, I * omega);

% compute omega_dot
omega_dot = I \ (tau_total - c);

omega_next = [];
q_dot = [];

% integrate omega if dt provided
if exist('dt','var') && ~isempty(dt)
    % simple RK4 integration for omega
    f = @(w, tau) I \ (tau - cross(w, I*w));
    k1 = f(omega, tau_total);
    k2 = f(omega + 0.5*dt*k1, tau_total);
    k3 = f(omega + 0.5*dt*k2, tau_total);
    k4 = f(omega + dt*k3, tau_total);
    omega_next = omega + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% compute quaternion derivative if requested
if integrate_q && ~isempty(q)
    % quaternion derivative q_dot (1x4) consistent with scalar-last convention
    % q_dot = 0.5 * Omega(omega) * q
    wx = omega(1); wy = omega(2); wz = omega(3);
    Omega = [  0   wz  -wy  wx;
              -wz   0   wx  wy;
               wy  -wx   0  wz;
              -wx  -wy  -wz  0 ];
    q_dot_col = 0.5 * Omega * q';  % 4x1
    q_dot = q_dot_col';
else
    q_dot = [];
end

end
