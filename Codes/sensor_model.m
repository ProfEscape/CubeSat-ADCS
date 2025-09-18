function [B_body_meas, sun_body_meas, gyro_meas] = sensor_model(q_body_from_eci, omega_body, r_eci_m, t_utc, params)
% SENSOR_MODEL - BMX160 magnetometer, gyroscope, and coarse sun sensors
%
% Inputs:
%   q_body_from_eci : quaternion (scalar-first, 4x1) Body <- ECI
%   omega_body      : true angular velocity [rad/s] in body frame
%   r_eci_m         : spacecraft position in ECI [m]
%   t_utc           : datetime UTC
%   params          : struct with sensor noise/std values from init_ADCS
%
% Outputs:
%   B_body_meas   : magnetometer reading [Tesla, body frame]
%   sun_body_meas : coarse sun sensor unit vector (body frame)
%   gyro_meas     : gyroscope reading [rad/s, body frame]

%% --- Magnetic field (IGRF) ---
% Convert ECI -> ECEF using GMST
jd = juliandate(t_utc);
gst = siderealTime(jd); % radians
R3 = [ cos(gst)  sin(gst) 0;
      -sin(gst)  cos(gst) 0;
       0         0        1 ];
r_ecef_m = R3 * r_eci_m;

% Geodetic latitude/longitude/altitude
[lat, lon, alt] = ecef2lla(r_ecef_m'); % Aerospace Toolbox, lat/lon in deg, alt in m

% Decimal year for IGRF
decYear = year(t_utc) + (day(t_utc,'dayofyear')-1 + hour(t_utc)/24 + ...
                         minute(t_utc)/1440 + second(t_utc)/86400)/365;

% Correct igrfmagm call: returns NED in nT
[XYZ, ~, ~, ~, ~] = igrfmagm(alt/1000, lat, lon, decYear);
B_ned_T = XYZ(:) * 1e-9; % Tesla

% Convert NED → ECEF


spheroid = wgs84Ellipsoid("kilometer"); % reference Earth model
[Bx_ecef, By_ecef, Bz_ecef] = ned2ecef( ...  %Mapping Toolbox
    B_ned_T(1), B_ned_T(2), B_ned_T(3), ...
    lat, lon, alt/1000, spheroid);

B_ecef = [Bx_ecef; By_ecef; Bz_ecef];

% Convert ECEF → ECI
B_eci = R3' * B_ecef';

% Rotate ECI → Body
R_body_from_eci = rotmat(q_body_from_eci, 'frame');
B_body_true = R_body_from_eci * B_eci;

% Add magnetometer noise (BMX160: 0.6 μT RMS)
B_body_meas = B_body_true + (0.6e-6)*randn(3,1);

%% --- Coarse sun sensor ---
% Sun vector in ECI (unit)
sun_eci = sunPosition(t_utc); % Aerospace Toolbox function
sun_body = R_body_from_eci * sun_eci;

% Six-panel photodiode model
normals = [ 1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1 ]'; % body axes
cosines = normals' * sun_body;
cosines(cosines<0) = 0;

% Albedo contribution (Earth Albedo toolbox)
albedo_percent = earth_albedo_model(r_ecef_m/1000, t_utc); % user toolbox  %%%%CHECK DOC
nadir_body = - R_body_from_eci * (r_eci_m / norm(r_eci_m));
albedo_cosines = max(0, normals' * nadir_body) * albedo_percent;

currents = cosines + albedo_cosines;
currents = currents + 0.01*randn(6,1); % shot noise

sun_vec_body = [ currents(1)-currents(2);
                 currents(3)-currents(4);
                 currents(5)-currents(6) ];

if norm(sun_vec_body) > 0
    sun_body_meas = sun_vec_body / norm(sun_vec_body);
else
    sun_body_meas = [0;0;0];
end

% Coarse sun sensor noise (±0.5° ≈ 8.7 mrad)
sun_body_meas = sun_body_meas + deg2rad(0.5)*randn(3,1);

%% --- Gyroscope (BMX160) ---
% True rate + bias + noise
% Datasheet: bias stability 3°/h = 0.000872 rad/s
% Output noise density 0.007 °/s/√Hz ≈ 1.22e-4 rad/s/√Hz
gyro_bias = (0.000872)*randn(3,1);     % random walk bias
gyro_noise = (1.22e-4/sqrt(params.sim.dt))*randn(3,1); % discrete-time noise
gyro_meas = omega_body + gyro_bias + gyro_noise;

end

%% --- Helpers ---
function gst = siderealTime(jd)
% Greenwich sidereal time (IAU 1982)
T = (jd - 2451545.0)/36525;
gmst = 280.46061837 + 360.98564736629*(jd - 2451545) ...
       + 0.000387933*T^2 - T^3/38710000;
gst = deg2rad(mod(gmst,360));
end

function sun_eci = sunPosition(t_utc)
% Sun position (unit vector, ECI) using Aerospace Toolbox
[rsun,~] = planetEphemeris(juliandate(t_utc),'Earth','Sun','421'); % km, ECI J2000
sun_eci = rsun'/norm(rsun);
end
