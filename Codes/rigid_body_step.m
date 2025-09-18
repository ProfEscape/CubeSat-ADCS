function [q_next, omega_next] = rigid_body_step(q, omega, torque, I, dt)
% Integrate attitude & angular velocity for rigid body using RK4
% q - quaternion scalar-first
% omega - angular velocity body frame (rad/s)
% torque - applied body torque [N*m]
% I - inertia matrix (3x3)
% dt - time step

% dynamics:
% omega_dot = I^{-1} * (torque - omega x (I*omega))
% qdot = 0.5 * Omega(omega) * q

% RK4 for omega and quaternion combined
state0.q = q; state0.omega = omega;

k1 = state_deriv(state0, torque, I);
k2 = state_deriv(state_add(state0, k1, dt/2), torque, I);
k3 = state_deriv(state_add(state0, k2, dt/2), torque, I);
k4 = state_deriv(state_add(state0, k3, dt), torque, I);

q_next = state0.q + (dt/6)*(k1.qdot + 2*k2.qdot + 2*k3.qdot + k4.qdot);
q_next = q_next / norm(q_next);
omega_next = state0.omega + (dt/6)*(k1.omega_dot + 2*k2.omega_dot + 2*k3.omega_dot + k4.omega_dot);

end

function s = state_add(s, k, factor)
s.q = s.q + factor * k.qdot;
s.omega = s.omega + factor * k.omega_dot;
end

function k = state_deriv(s, torque, I)
omega = s.omega;
Iomega = I * omega;
omega_dot = I \ (torque - cross(omega, Iomega));
Omega = [0    -omega(1) -omega(2) -omega(3);
         omega(1) 0     omega(3) -omega(2);
         omega(2) -omega(3) 0     omega(1);
         omega(3) omega(2) -omega(1) 0];
qdot = 0.5 * Omega * s.q;
k.qdot = qdot;
k.omega_dot = omega_dot;
end
