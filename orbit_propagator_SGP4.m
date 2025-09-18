function [r_eci_km, v_eci_kms] = orbit_propagator_SGP4(tle_line1, tle_line2, epoch)
% ORBIT_PROPAGATOR_SGP4  Propagate TLE to ECI r,v using Vallado's SGP4 MATLAB port
%
% Inputs:
%   tle_line1, tle_line2 - strings, the TLE lines
%   epoch                - datetime (UTC) or Julian date (numeric)
%
% Outputs:
%   r_eci_km  - 3x1 ECI position [km]
%   v_eci_kms - 3x1 ECI velocity [km/s]

whichconst = 'wgs72';   % Vallado recommended for TLEs
typerun    = 'c';       % catalog mode
typeinput  = 'e';       % use epoch input

% Initialize satellite record
[satrec, ~, ~, ~] = twoline2rv(whichconst, tle_line1, tle_line2, typerun, typeinput);

% Determine current epoch in Julian date
if isa(epoch,'datetime')
    jd_epoch = juliandate(epoch);
else
    jd_epoch = epoch;
end

% Compute tsince [min] since TLE epoch
tsince_min = (jd_epoch - satrec.jdsatepoch) * 24 * 60;

% Propagate with SGP4
[~, r_eci_km, v_eci_kms] = sgp4(satrec, tsince_min);

% Ensure column vectors
r_eci_km  = r_eci_km(:);
v_eci_kms = v_eci_kms(:);
end
