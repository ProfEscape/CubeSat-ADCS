%%Helper Function (MEKF)
function R = quat_to_rotm(q)
% QUAT_TO_ROTM  Convert a scalar-last quaternion q = [qx qy qz qw]
% to the rotation matrix R that maps vectors from inertial -> body:
% r_body = R * r_inertial
%
% Note: this assumes the same convention used in the rest of the code:
% quaternion is 1x4 with scalar-last: [qx qy qz qw]

% ensure normalized
q = quatnormalize(q);

qx = q(1); qy = q(2); qz = q(3); qw = q(4);

% rotation matrix (scalar-last convention)
R = [ 1-2*(qy^2+qz^2),   2*(qx*qy - qz*qw),  2*(qx*qz + qy*qw);
      2*(qx*qy + qz*qw), 1-2*(qx^2+qz^2),    2*(qy*qz - qx*qw);
      2*(qx*qz - qy*qw), 2*(qy*qz + qx*qw),  1-2*(qx^2+qy^2) ];
end
