function [q_pred, b_pred, P_pred, Phi_k, Qs] = mekf_predict(q_hat, b_hat, P, w_meas, dt, Qproc)
% MEKF_PREDICT - Multiplicative EKF predict step (DDJF Algorithm 4)
%
% Inputs:
%   q_hat  - 1x4 estimated quaternion [qx qy qz qw] at time k
%   b_hat  - 3x1 estimated gyro bias
%   P      - 6x6 covariance at time k
%   w_meas - 3x1 gyro measurement (rad/s)
%   dt     - timestep (s)
%   Qproc  - 6x6 process noise covariance
%
% Outputs:
%   q_pred - 1x4 predicted quaternion
%   b_pred - 3x1 predicted bias
%   P_pred - 6x6 predicted covariance
%   Phi_k  - 6x6 discrete state transition matrix
%   Qs     - 6x6 discrete process noise contribution

% Bias-corrected angular velocity
w_hat = w_meas - b_hat;

% Quaternion derivative (Omega matrix form)
Omega = [  0      -w_hat(1) -w_hat(2) -w_hat(3);
           w_hat(1)     0    w_hat(3) -w_hat(2);
           w_hat(2) -w_hat(3)    0     w_hat(1);
           w_hat(3)  w_hat(2) -w_hat(1)    0  ];

q_dot = 0.5 * (Omega * q_hat')';
q_pred = q_hat + q_dot * dt;
q_pred = quatnormalize(q_pred);

% Bias unchanged
b_pred = b_hat;

% Continuous-time system matrix
F = [ -skew3(w_hat), -eye(3);
       zeros(3),      zeros(3) ];

% Block matrix exponential
A = [-F,          Qproc;
      zeros(6),   F'];
B = expm(A*dt);

Phi_k = B(7:12,7:12);
Qs    = Phi_k * B(1:6,7:12);

% Covariance propagation
P_pred = Phi_k * P * Phi_k' + Qs;
end
