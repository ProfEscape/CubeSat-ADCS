function sun_eci = sun_vector_vallado(t_utc)
% sun_vector_vallado - compute Sun unit vector in ECI at given UTC datetime.
% Implemented with the standard (Vallado-like) low-cost ephemeris (sufficient for ADCS).
% Input: t_utc - MATLAB datetime (UTC)
% Output: sun_eci - 3x1 unit vector (ECI frame, same convention used elsewhere)

% Convert to Julian centuries since J2000
jd = juliandate(t_utc);
T = (jd - 2451545.0)/36525;

% Mean anomaly (deg)
M = 357.5277233 + 35999.05034 * T;
% Mean longitude (deg)
L = 280.460 + 36000.771 * T;
% Eccentricity of Earth's orbit
e = 0.016708617 - 0.000042037 * T - 0.0000001236 * T^2;

% Sun's equation of center
M_rad = deg2rad(mod(M,360));
C = (1.914602 - 0.004817*T - 0.000014*T^2)*sin(M_rad) + ...
    (0.019993 - 0.000101*T)*sin(2*M_rad) + 0.000289*sin(3*M_rad);

% True longitude (deg)
lambda = L + C;
lambda_rad = deg2rad(mod(lambda,360));

% Obliquity of ecliptic (deg)
eps = 23.439291 - 0.0130042*T - 1.64e-7*T^2 + 5.04e-7*T^3;
eps_rad = deg2rad(eps);

% Sun vector in ecliptic coordinates (assuming unit distance)
x = cos(lambda_rad);
y = sin(lambda_rad);
z = 0;

% Rotate to ECI
R = [1, 0, 0; 0, cos(eps_rad), -sin(eps_rad); 0, sin(eps_rad), cos(eps_rad)];
sun_eci = R * [x; y; z];
sun_eci = sun_eci / norm(sun_eci);
end
