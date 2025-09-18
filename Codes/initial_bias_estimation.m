function b0 = initial_bias_estimation(q_seq, gyro_seq, dt)
% INITIAL_BIAS_ESTIMATION - Estimate initial gyro bias using Wahba quaternions
%
% Inputs:
%   q_seq    - n x 4 quaternion array [qx qy qz qw], Wahba solutions
%   gyro_seq - 3 x n gyro measurements
%   dt       - timestep (s)
%
% Output:
%   b0       - 3x1 estimated gyro bias

n = size(q_seq,1);
if n < 2
    error('Need at least 2 quaternions');
end

omega_rel = zeros(3,n-1);
for k = 2:n
    dq = quatmultiply(q_seq(k,:), quatinv(q_seq(k-1,:)));
    dq = quatnormalize(dq);

    % Convert quaternion -> axis-angle -> rotation vector
    axang = quat2axang(dq); % [ax ay az angle]
    phi = axang(1:3) * axang(4); % rotation vector
    omega_rel(:,k-1) = phi(:)/dt;
end

omega_mean = mean(omega_rel,2);
gyro_mean  = mean(gyro_seq(:,2:n),2);

b0 = omega_mean - gyro_mean;
end
