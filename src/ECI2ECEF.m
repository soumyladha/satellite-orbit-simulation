%*************************************************************************
% This function gives the coordinate transformation matrix from J2000 ECI 
% to ECEF 
% Function Argument :
%                     date
%                     time 
% Library calls : NIL

% Function Outputs : 
%                   c_eci_to_ecef
% Functions calls : 
%                  nutation
%                  precession
%                  rotz
%                  GMST
% Global Variables : NIL
% % Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************
function c_eci_to_ecef = ECI2ECEF(date_time,deltaUT1,x_pole,y_pole,tt_utc)

date_d = date_time+(tt_utc)/86400;

precession_matrix = Precession(date_d);
precession_matrix = precession_matrix';

[nutation_matrix,psi,e] = Nutation(date_d);

% Earth rotation For reference refer Satellite Orbits, Pg 181, Eqn 5.66
theta = GMST(date_time,deltaUT1)+ psi*cos(e);

earth_rotation_matrix = Rotz(theta);

% Conversion from seconds to radians

x_pole=x_pole*pi/(3600*180);
y_pole=y_pole*pi/(3600*180);
polar_motion_matrix = [1 0 x_pole;0 1 -y_pole;-x_pole y_pole 1];

c_eci_to_ecef = polar_motion_matrix*earth_rotation_matrix*nutation_matrix*precession_matrix;


