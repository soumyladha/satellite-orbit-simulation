%*************************************************************************
% This function gives GMST for julian date 
% Function Argument :
%                     Julian date , delta UT1 
% Library calls : NIL

% Function Outputs : 
%                   GMST
% Functions calls : 
%                  
% Global Variables : NIL
% % Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************
function gree_side_time = GMST(julian_date,delta_ut1)

store_date = julian_date;
julian_date = julian_date+delta_ut1;
% Calculating julian date at midnight of given date
julian_date_min = floor(julian_date)-0.5;
julian_date_max = floor(julian_date)+0.5;

if julian_date_max<=julian_date
    julian_date_UT0 = julian_date_max;
elseif julian_date>julian_date_min && julian_date<julian_date_max
    julian_date_UT0 = julian_date_min;   
end
julian_date = store_date;
juliant_0 = (julian_date_UT0-2451545)/36525;
juliant = (julian_date-2451545)/36525;
juliant_1 = (juliant-juliant_0)*36525*86400+delta_ut1;

% GMST converts seconds,minutes and hours to hours 
% For reference please refer Satellite Orbits, Pg 167, Eqn 5,19
gree_side_time = 24110.54841+8640184.812866*juliant_0...
    +1.002737909350795*juliant_1+0.093104*juliant*juliant-...
    0.0000062*juliant*juliant*juliant;

% GMST in radians
gree_side_time= 2*pi*mod(gree_side_time,86400)/86400;


