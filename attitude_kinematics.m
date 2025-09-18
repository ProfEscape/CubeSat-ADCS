function q_next = attitude_kinematics(q, omega_body, dt)
% ATTITUDE_KINEMATICS  Propagate quaternion using body angular velocity
%
% Inputs:
%   q           - 1x4 quaternion [qx qy qz qw] (scalar-last)
%   omega_body  - 3x1 angular velocity in body frame (rad/s)
%   dt          - time-step (s)
%
% Output:
%   q_next      - 1x4 quaternion after dt (normalized)
%
% Kinematics used (DDJF 6.2.4.2 style):
%   q_dot = 0.5 * Omega(omega) * q
% with Omega built for scalar-last convention.
%
% Integration: RK4 for accuracy.

% helper - quaternion derivative
    function qd = qdot(q_in, w)
        % q_in as 1x4 [qx qy qz qw], w as 3x1
        % Build Omega in scalar-last frame (so that q_dot = 0.5 * Omega * q')
        wx = w(1); wy = w(2); wz = w(3);
        % Omega acting on column quaternion [qx qy qz qw]'
        Omega = [  0   wz  -wy  wx;
                  -wz   0   wx  wy;
                   wy  -wx   0  wz;
                  -wx  -wy  -wz  0 ];
        qd_col = 0.5 * Omega * q_in';  % 4x1
        qd = qd_col';
    end

% RK4 integration
k1 = qdot(q, omega_body);
q1 = q + 0.5 * dt * k1;
k2 = qdot(q1, omega_body);
q2 = q + 0.5 * dt * k2;
k3 = qdot(q2, omega_body);
q3 = q + dt * k3;
k4 = qdot(q3, omega_body);

q_next = q + (dt/6) * (k1 + 2*k2 + 2*k3 + k4);

% normalize using Aerospace toolbox if present, else manual
if exist('quatnormalize','file') == 2
    q_next = quatnormalize(q_next);
else
    q_next = q_next / norm(q_next);
end
end
