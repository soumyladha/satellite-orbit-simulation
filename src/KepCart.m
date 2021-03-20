%*************************************************************************
% This function calculates the pos, vel ,acc and jerk of a satellite
%  from Kepler's orbital parameters for a given time instant
% 
% Function Calls :  Rotx (Rotation Matrix about X)
%                   Rotz (Rotation Matrix about Z)
%                   EccenAnom (Calculation of Eccentric from Mean Anomaly)
% Function Arguments:
%                   semi_axis = Semi major axis
%                   eccen = Eccentricty of Orbit
%                   incli = Inclination
%                   acend_node = The Right ascension of the ascending node in degree
%                   arg_per = argument of perigee in degree
%                   mean_anom = Mean anomaly in degree
%                   flag = whether the inputs are in Inertial or ECEF frame
%                   JulianDate = value of date in Julian
%                    
% Return Value :
%                   
%                    posPar = position
%                    velPar = velocity
%                    accPar = accelaration
%
% Library Calls :
%                   cross (Taking cross product of two vectors)
%                   norm (Taking magnitude of vector)
%                   dot (Taking dot product of two vectors)
%                   atan (Taking tan inverse)
%                   sqrt (square root) 
% Global Variables : NIL
% 
% 
% Version History: 
%                   <1.1> <Soumy Ladha>
%**************************************************************************

function [pos_par,vel_par,acc_par,jerk_par] = KepCart(semi_axis,eccen,incli,...
    acend_node,arg_per,mean_anom,flag,~)

if flag == 1% orbital elements in the J2000 ECI frame
    Constants; % Loading Constants from file

    % Calculation of Mean Motion
    % For reference refer book Satellite Orbits Page 23 ,Equation 2.35
    mean_m = sqrt(Grav_C*M_earth/semi_axis^3); 

    % Calculation of Eccentric Anomaly
    % For reference refer book Satellite Orbits Page 24 ,Equation 2.42
    ecc_anom = EccenAnom(eccen,mean_anom);

    % Calculation of Rate of change of eccentric anomaly
    % For reference refer book Satellite Orbits Page 23 ,Equation 2.34
    ecc_anom_dot = mean_m/(1-eccen*cos(ecc_anom));

    % Calculation of position vector in Orbital plane
    % For reference refer book Satellite Orbits Page 22 ,Equation 2.30
    pos_par_orb_x = semi_axis*(cos(ecc_anom)-eccen); 
    pos_par_orb_y = semi_axis*sqrt(1-eccen^2)*sin(ecc_anom);
    pos_par_orb_z = 0;

    % Combining into a vector
    pos_par_orb = [pos_par_orb_x; pos_par_orb_y; pos_par_orb_z];

    % Calculation of velocity Vector in Orbital Plane
    % For reference refer attached Document (will be attached later)
    vel_par_orb_x = semi_axis*ecc_anom_dot*-sin(ecc_anom);
    vel_par_orb_y = semi_axis*ecc_anom_dot*sqrt(1-eccen^2)*cos(ecc_anom);
    vel_par_orb_z = 0;
    
    % Combining into a vector
    vel_par_orb = [vel_par_orb_x;vel_par_orb_y;vel_par_orb_z];

    % Calculation of accelaration in Orbital Plane
    % For reference refer attached Document (will be attached later)
    accpar_orb_x = (ecc_anom_dot^2)*semi_axis*(eccen-cos(ecc_anom))/(1-eccen*cos(ecc_anom));
    accpar_orb_y = (ecc_anom_dot^2)*semi_axis*sqrt(1-eccen^2)*(-sin(ecc_anom))/(1-eccen*cos(ecc_anom));
    accpar_orb_z = 0;

    % Combining into a vector
    accpar_orb = [accpar_orb_x;accpar_orb_y;accpar_orb_z];

    % Calculation of jerk in Orbital Plane
    % For reference refer attached Document (will be attached later)
    % Defining new variable
    grav_const_mass_earth = Grav_C*M_earth;
     pos_par_orb_mag = sqrt(sum(pos_par_orb.*pos_par_orb));
     pos_par_orb_mag_dot = -1*semi_axis*eccen*sin(ecc_anom)*ecc_anom_dot;
     
    jerkpar_orb_x =-3*grav_const_mass_earth*pos_par_orb_mag_dot*pos_par_orb_x/(pos_par_orb_mag^4)...
        -grav_const_mass_earth*vel_par_orb_x/(pos_par_orb_mag^3);
    
    jerkpar_orb_y=-3*grav_const_mass_earth*pos_par_orb_mag_dot*pos_par_orb_y/(pos_par_orb_mag^4)...
        -grav_const_mass_earth*vel_par_orb_y/(pos_par_orb_mag^3);
    
    jerkpar_orb_z=-3*grav_const_mass_earth*pos_par_orb_mag_dot*pos_par_orb_z/(pos_par_orb_mag^4)...
        -grav_const_mass_earth*vel_par_orb_z/(pos_par_orb_mag^3);

    % Combining into a vector
    jerkpar_orb = [jerkpar_orb_x;jerkpar_orb_y;jerkpar_orb_z];

    % Rotation matrix form transformation from Orbital to Cartesian
    % For reference refer to the book Satellite Orbits Page 27 ,Equation
    % 2.5
    rotation = Rotz(-acend_node)*Rotx(-incli)*Rotz(-arg_per);

    % position vector in Cartesian
    pos_par = (rotation*pos_par_orb);

    % velocity vector in Cartesian
    vel_par = (rotation*vel_par_orb);

    % accelaration vector in Cartesian
    acc_par = (rotation*accpar_orb);

    % jerk vector in Cartesian
    jerk_par = (rotation*jerkpar_orb);
    
    
    pos_par = pos_par';
    vel_par = vel_par';
    acc_par = acc_par';
    jerk_par = jerk_par';
    
    
else % Orbital elements in the ECEF frame at the input time
    Constants; % Loading Constants from file

    % Calculation of Mean Motion
    % For reference refer book Satellite Orbits Page 23 ,Equation 2.35
    mean_m = sqrt(Grav_C*M_earth/semi_axis^3); 

    % Calculation of Eccentric Anomaly
    % For reference refer book Satellite Orbits Page 24 ,Equation 2.42
    ecc_anom = EccenAnom(eccen,mean_anom);

    % Calculation of Rate of change of eccentric anomaly
    % For reference refer book Satellite Orbits Page 23 ,Equation 2.34
    ecc_anom_dot  =  mean_m/(1-eccen*cos(ecc_anom));
    
    % Calculation of position vector in Orbital plane
    % For reference refer book Satellite Orbits Page 22 ,Equation 2.30
    pos_par_orb_x = semi_axis*(cos(ecc_anom)-eccen); 
    pos_par_orb_y = semi_axis*sqrt(1-eccen^2)*sin(ecc_anom);
    pos_par_orb_z = 0;
    
    % Combining into the vector
    pos_par_orb = [pos_par_orb_x; pos_par_orb_y; pos_par_orb_z];
    
    % pos_par_orb_mag calculates the magnitude of the position
    pos_par_orb_mag = sqrt(sum(pos_par_orb.*pos_par_orb));
    
    % pos_par_orb_mag_dot calculates the rate of change of magnitude
    % of position
    pos_par_orb_mag_dot = -1*semi_axis*eccen*sin(ecc_anom)*ecc_anom_dot;
    
    % Calculation of velocity Vector in Orbital Plane
    % For reference refer attached Document (will be attached later)
    vel_par_orb_x = semi_axis*ecc_anom_dot*-sin(ecc_anom);
    vel_par_orb_y = semi_axis*ecc_anom_dot*sqrt(1-eccen^2)*cos(ecc_anom);
    vel_par_orb_z = 0;

    % Combining into the vector
    vel_par_orb = [vel_par_orb_x;vel_par_orb_y;vel_par_orb_z];
  
    % GravConsXMassEarth is the product of Universal Gravitational Constant
    % & Mass of earth
    grav_const_mass_earth = Grav_C*M_earth;
    
    % Calculation of accelaration in Orbital Plane
    % For reference refer attached Document (will be attached later)
    accpar_orb_x = -grav_const_mass_earth*pos_par_orb_x/(pos_par_orb_mag^3);
    accpar_orb_y = -grav_const_mass_earth*pos_par_orb_y/(pos_par_orb_mag^3);
    accpar_orb_z = -grav_const_mass_earth*pos_par_orb_z/(pos_par_orb_mag^3);
    
    accpar_orb = [accpar_orb_x;accpar_orb_y;accpar_orb_z];

    
    % Calculation of jerk in Orbital Plane
    % For reference refer attached Document (will be attached later)
    jerkpar_orb_x =-3*grav_const_mass_earth*pos_par_orb_mag_dot*pos_par_orb_x/(pos_par_orb_mag^4)...
        -grav_const_mass_earth*vel_par_orb_x/(pos_par_orb_mag^3);
    
    jerkpar_orb_y=-3*grav_const_mass_earth*pos_par_orb_mag_dot*pos_par_orb_y/(pos_par_orb_mag^4)...
        -grav_const_mass_earth*vel_par_orb_y/(pos_par_orb_mag^3);
    
    jerkpar_orb_z=-3*grav_const_mass_earth*pos_par_orb_mag_dot*pos_par_orb_z/(pos_par_orb_mag^4)...
        -grav_const_mass_earth*vel_par_orb_z/(pos_par_orb_mag^3);
    
    jerkpar_orb = [jerkpar_orb_x;jerkpar_orb_y;jerkpar_orb_z];
    
    % Rotation matrix form transformation Orbital to Cartesian
    % For reference refer to the book Satellite Orbits Page 27 ,Equation
    % 2.5
    rotation = Rotz(-acend_node)*Rotx(-incli)*Rotz(-arg_per);
    
    % position vector in Cartesian frame
    pos_cart = (rotation*pos_par_orb);

    % velocity vector in Cartesian frame
    vel_cart = (rotation*vel_par_orb);

    % accelaration vector in Cartesian frame
    acc_cart = (rotation*accpar_orb);

    % jerk vector in Cartesian Inertial
    jerk_cart = (rotation*jerkpar_orb);
   
    
    % Additional terms for ECEF frame
    
    % position vector in ECEF
    % For reference see the attached document
    pos_par = pos_cart;
    
    % velocity vector in ECEF
    % For reference see the attached document
    vel_par = vel_cart - cross(Omega_earth,pos_cart);
    
    % accelaration vector in ECEF
    % For reference see the attached document
    acc_par = acc_cart - 2*cross(Omega_earth,vel_par)...
                        - cross(Omega_earth,cross(Omega_earth,pos_cart));
                    
    % Jerk vector in ECEF
    % For reference see the attached document
    jerk_par = jerk_cart - 3*cross(Omega_earth,acc_par)...
                - 3*cross(Omega_earth,cross(Omega_earth,vel_par))...
     - cross(Omega_earth,cross(Omega_earth,cross(Omega_earth,pos_cart)));
    
    pos_par = pos_par';
    vel_par = vel_par';
    acc_par = acc_par';
    jerk_par = jerk_par';
    
    
    
end

