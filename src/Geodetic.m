%*************************************************************************
% The following function calculates the Geodetic latitude, longitude and
% height from postion vector
% Function Arguments : x_pos
%                      y_pos
%                      z_pos
% Functions Outputs : Latitude (degrees)
%                     Longitude (degrees)
%                     height (m)
% Library calls : NIL
% Global Variables : NIL
% Version History:
%                   <1.1> <Soumy Ladha>
%*************************************************************************
function [latitude,longitude,height] = Geodetic(x_pos,y_pos,z_pos)
temp = sqrt(x_pos^2+y_pos^2);
if temp < 10^(-7) % For polar coordinate correction
    if z_pos <0
        latitude = -90;
        longitude = 0;
        height = z_pos+semi_minor;
    else
        latitude = 90;
        longitude = 0;
        height = z_pos-semi_minor;
    end   
else
    
    longitude = atan2(y_pos,x_pos); % Initial Guess
    latitude = atan(z_pos/temp); % Initial Guess
    
    Constants;
    
    eccentricity = sqrt(2*flatness-flatness^2); % For refrence GNSS,inertial
    % and multisensor navigation by Paul D Groves, Page 37, Equation 2.53
    
    height = 0; % Initialization
    height_old = sqrt(x_pos*x_pos+y_pos*y_pos+z_pos*z_pos)-6378137; % Initialization
    
    % For reference, refer GNSS,inertial and multisensor navigation by Paul
    % D Groves, Page 38
    while abs(height_old-height)>10^-4 % Stopping condition Error tolerance
        r_n_phi = (semi_major/sqrt(1-eccentricity^2*(sin(latitude))^2));
        latitude = atan((z_pos/temp)*(1-(eccentricity^2)*(r_n_phi)/(r_n_phi+height))^(-1));
        height_old = height;
        height = (temp/cos(latitude))-r_n_phi;    
    end   
    
    % Converting to radians
    latitude = latitude*180/pi;
    longitude = longitude*180/pi;
end


