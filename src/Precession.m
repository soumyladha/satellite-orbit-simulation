%*************************************************************************
% This function gives the precission matrix for the given Julian date 
% Function Argument :
%                     julian_date 
% Library calls : NIL

% Function Outputs : 
%                   precession
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

function precession_matrix = Precession(date)
    
    
    
    % Time in julian centuries
    t_j_current = (date-2451545.0)/36525;
    
    % For reference please refer Precession matrix absed on IAU (1976) System of Astronomical Constants
    t_j2000 = 0; % For specific case when epoch 1 is January 1 2000
    % Coordinate change due to Precession
    eta_a = ((DegMinSec2Degree(0,0,2306.2181) + DegMinSec2Degree(0,0,1.39656)*t_j2000...
        -DegMinSec2Degree(0,0,0.000139)*t_j2000*t_j2000))*t_j_current...
        +(DegMinSec2Degree(0,0,0.30188)-DegMinSec2Degree(0,0,0.000344)*t_j2000)*t_j_current*t_j_current...
        + DegMinSec2Degree(0,0,0.017998)*t_j_current*t_j_current*t_j_current;
    
    zeta_a = ((DegMinSec2Degree(0,0,2306.2181) + DegMinSec2Degree(0,0,1.39656)*t_j2000...
        -DegMinSec2Degree(0,0,0.000139)*t_j2000*t_j2000))*t_j_current...
        +(DegMinSec2Degree(0,0,1.09468)+DegMinSec2Degree(0,0,0.000066)*t_j2000)*t_j_current*t_j_current...
        + DegMinSec2Degree(0,0,0.018203)*t_j_current*t_j_current*t_j_current;
    
    theta_a = ((DegMinSec2Degree(0,0,2004.3109) - DegMinSec2Degree(0,0,0.85330)*t_j2000...
        -DegMinSec2Degree(0,0,0.000217)*t_j2000*t_j2000))*t_j_current...
        -(DegMinSec2Degree(0,0,0.42665)-DegMinSec2Degree(0,0,0.000217)*t_j2000)*t_j_current*t_j_current...
        - DegMinSec2Degree(0,0,0.041833)*t_j_current*t_j_current*t_j_current;
    
    eta_a = eta_a*pi/180;
    zeta_a = zeta_a*pi/180;
    theta_a = theta_a*pi/180;
    
    precession_matrix = Rotz(-zeta_a)*Roty(theta_a)*Rotz(-eta_a);
    precession_matrix = precession_matrix';
     

end