function [q_upd, b_upd, P_upd, eps_k] = mekf_correct(q_pred, b_pred, P_pred, meas, refs, R)
% MEKF_CORRECT - MEKF correction step per DDJF
%
% Inputs:
%   q_pred  - 1x4 predicted quaternion [qx qy qz qw] (scalar-last)
%   b_pred  - 3x1 predicted gyro bias
%   P_pred  - 6x6 predicted covariance
%   meas    - cell array {z1, z2, ...} body-frame measurements (3x1 each)
%   refs    - cell array {r1, r2, ...} inertial reference vectors (3x1 each)
%   R       - measurement covariance (3n x 3n)
%
% Outputs:
%   q_upd   - updated quaternion (1x4)
%   b_upd   - updated bias (3x1)
%   P_upd   - updated covariance (6x6)
%   eps_k   - error-state correction vector [dtheta; dbias] (6x1)

% Number of vector measurements and sizes
numVec = numel(meas);
m = 3 * numVec;

% Allocate stacked measurement / prediction / Jacobian
z = zeros(m,1);
h = zeros(m,1);
H = zeros(m,6);

% Compute rotation matrix A(q_pred) that maps inertial -> body
A_q = quat_to_rotm(q_pred);   % 3x3

for k = 1:numVec
    idx = (k-1)*3 + (1:3);
    r_eci = refs{k}(:);       % 3x1 inertial/reference vector

    % predicted measurement: r_body = A(q_pred) * r_eci
    b_hat_meas = (A_q * r_eci);   % 3x1

    % stack predicted measurement and actual measurement
    h(idx) = b_hat_meas;
    z(idx) = meas{k}(:);

    % EXACT analytic Jacobian per DDJF:
    % H_i = [ A(q_hat) * skew(r_i) ,  0_3x3 ]
    H(idx,1:3) = A_q * skew3(r_eci);
    H(idx,4:6) = zeros(3,3);
end

% Kalman gain
S = H * P_pred * H' + R;
K = (P_pred * H') / S;  % 6 x m

% innovation
nu = z - h;

% error-state correction
eps_k = K * nu;    % 6x1: [dtheta; dbias]
dtheta = eps_k(1:3);
dbias  = eps_k(4:6);

% multiplicative quaternion update (small-angle)
% represent small rotation as quat dq = [0.5*dtheta, 1] (scalar-last)
dq = [ (0.5*dtheta(:))'  1 ];   % 1x4
dq = quatnormalize(dq);

% update quaternion: q_new = dq ⊗ q_pred
q_upd = quatmultiply(dq, q_pred);
q_upd = quatnormalize(q_upd);

% update bias
b_upd = b_pred + dbias;

% covariance update (Joseph form for numerical stability)
I6 = eye(6);
P_upd = (I6 - K*H) * P_pred * (I6 - K*H)' + K * R * K';

end
