%*************************************************************************
% This function gives the Nutation matrix for the given Julian date 
% Function Argument :
%                     date 
% Library calls : NIL

% Function Outputs : 
%                   nutation,
%                   delta_final_psi
%                   Angle between Equatorial and Ecliptic Plane = angle_eq_ecl
% Functions calls : 
%                   DegMinSec2Degree -> Converts Degrees, Minutes and Seconds to
%                   Degrees
%                   Rotz
%                   Roty
%                   Rotx
% Global Variables : NIL
% % Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************
function [nutation_matrix,delta_final_psi,angle_eq_ecl] = Nutation(date)

% Time in julian centuries
time = (date-2451545.0)/36525; 

% For reference refer Satellite orbits Pg 178 5.3.3
moon_anomaly = DegMinSec2Degree(134,57,46.733)+DegMinSec2Degree(477198,52,02.633)*time...
    +DegMinSec2Degree(0,0,31.310)*time*time+DegMinSec2Degree(0,0,0.064)*time*time*time;

sun_anomaly = DegMinSec2Degree(357,31,39.804)+DegMinSec2Degree(35999,03,01.244)*time...
    -DegMinSec2Degree(0,0,0.577)*time*time-DegMinSec2Degree(0,0,0.012)*time*time*time;

mean_moon_distance = DegMinSec2Degree(93,16,18.877)+DegMinSec2Degree(483202,01,03.137)*time...
    -DegMinSec2Degree(0,0,13.257)*time*time+DegMinSec2Degree(0,0,0.011)*time*time*time;

delta_mean_longitudes = DegMinSec2Degree(297,51,01.307)+DegMinSec2Degree(445267,06,41.328)*time...
    -DegMinSec2Degree(0,0,6.891)*time*time+DegMinSec2Degree(0,0,0.019)*time*time*time;

acc_node_lunar= DegMinSec2Degree(125,02,40.280)-DegMinSec2Degree(1934,08,10.539)*time...
    +DegMinSec2Degree(0,0,7.455)*time*time+DegMinSec2Degree(0,0,0.008)*time*time*time;

nut80=[];
Constants % importing data of nutation table

phi_constants = [moon_anomaly sun_anomaly mean_moon_distance...
    delta_mean_longitudes acc_node_lunar];
phi_constants = phi_constants*pi/180;

phi_table = nut80(:,1:5);
phi_vector = phi_table*phi_constants';

delta_psi = nut80(:,6:7)*[1; time];
delta_e = nut80(:,8:9)*[1; time];

delta_psi= delta_psi';
delta_e = delta_e';

% For Reference please refer Satellite Orbits, Pg 178, Equation 5.58
delta_final_psi = delta_psi*sin(phi_vector)*(0.0001/3600)*pi/180;
delta_final_e = delta_e*cos(phi_vector)*(0.0001/3600)*pi/180;

% For reference refer satellite orbits equation 5.42
angle_eq_ecl=DegMinSec2Degree(23.4392911,0,0)-DegMinSec2Degree(0,0,46.8150)*time...
    -DegMinSec2Degree(0,0,0.00059)*time*time+DegMinSec2Degree(0,0,0.001813)*time*time*time;

angle_eq_ecl=angle_eq_ecl*pi/180;
nutation_matrix =  Rotx(-angle_eq_ecl-delta_final_e)*Rotz(-delta_final_psi)*Rotx(angle_eq_ecl);



