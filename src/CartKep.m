%*************************************************************************
% The following function calculates the 6 kepler's orbital parameters from
% position and velocity given in cartesian coordinates
% Functions Called :  
%                   NIL
% Library Calls :
%                   Cross (Taking cross product of two vectors)
%                   norm (Taking magnitude of vector)
%                   dot (Taking dot product of two vectors)
%                   atan (Taking tan inverse)
%                   sqrt (square root) 
% Functions input : 
%                   Position in m  
%                   Velocity in m/s
%                   Flag if the input is in ECI , ECEF
%                   julian_date = value of date in Julian
% Function Output :
%                   semi_axis = Semi major axis
%                   eccen = Eccentricty of Orbit
%                   incli = Inclination
%                   acend_node = The Right ascension of the ascending node 
%                                in degrees
%                   arg_per = argument of perigee in degrees
%                   mean_anom = Mean anomaly in degrees
% Global Variables : NIL
% 
% 
% Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************
function [semi_axis,eccen,incli,acend_node,arg_per,mean_anom]=CartKep(position,...
    velocity,flag,~)


if flag == 1 % Calculation in ECI frame
    
     Constants; %loading constants from file
     
     

     % Calculation of angular momentum
     ang_mom = cross(position,velocity);
     
     
        
    % Calculation of inclination
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.58
    if ang_mom(3)>0
        incli=atan(sqrt(ang_mom(1)^2+ang_mom(2)^2)/ang_mom(3));
    else
        incli=pi+atan(sqrt(ang_mom(1)^2+ang_mom(2)^2)/ang_mom(3));
    end
    incli=incli*180/pi;

    % Calculation of The Right ascension of the ascending node
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.58
    if ang_mom(2)<0
        if(ang_mom(1)>0)
            acend_node=atan(ang_mom(1)/-ang_mom(2));
        else
            acend_node=2*pi+atan(ang_mom(1)/-ang_mom(2));
            
        end
    else
        acend_node=pi+atan(ang_mom(1)/-ang_mom(2));
    end
    acend_node=acend_node*180/pi;

    % Calculation of semi latus rectum
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.59
    sem_latus = (norm(ang_mom))^2/(Grav_C*M_earth);

    % Calculation of Semi major axis
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.60
    semi_axis = 1/(2/norm(position)-(norm(velocity)^2)/(Grav_C*M_earth));

    % Calculatio of Mean motion
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.61
    mean_m = sqrt(Grav_C*M_earth/semi_axis^3); 

    % Calculation of Eccentricity
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.62
    eccen = sqrt(1-(sem_latus/semi_axis));

    % Calculation of Eccentric Anomaly
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.64
    if (1-norm(position)/semi_axis)>0
        ecc_anom = atan(dot(position,velocity)/(semi_axis^2*mean_m*(1-norm(position)/semi_axis)));
    else
        ecc_anom = pi+ atan(dot(position,velocity)/(semi_axis^2*mean_m*(1-norm(position)/semi_axis)));
    end

    % Calculation of Mean Anomaly
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.65
	mean_anom = ecc_anom-eccen*sin(ecc_anom);
    if mean_anom>=0
        mean_anom=mean_anom*180/pi;
    else
        mean_anom=mean_anom*180/pi+360;        
    end
    
    % Calculation of Argument of perigee
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.66
    % Calculation of Unit vector Along Anglar momemtum Vector
    ang_mom_u=ang_mom/norm(ang_mom);
    % Calculation of Argument of Latitude
    if (-position(1)*ang_mom_u(2)+position(2)*ang_mom_u(1))>0
        arg_lat = atan(position(3)/(-position(1)*ang_mom_u(2)+position(2)*ang_mom_u(1)));
    else
        arg_lat =pi+ atan(position(3)/(-position(1)*ang_mom_u(2)+position(2)*ang_mom_u(1)));
    end

    % Calculation of True Anomaly
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.67
    if (cos(ecc_anom)-eccen)>0
        true_anom = atan(sqrt(1-eccen^2)*sin(ecc_anom)/(cos(ecc_anom)-eccen));
    else
        true_anom=pi+atan(sqrt(1-eccen^2)*sin(ecc_anom)/(cos(ecc_anom)-eccen));
    end

    arg_per = arg_lat-true_anom;
    if (arg_per>0)
        arg_per = arg_per*180/pi;  
    else
        arg_per = arg_per*180/pi+360;
    end
else
    
    Constants; %loading constants from file
    
    velocity = velocity + cross(Omega_earth,position);%calculation of velocity in inertial frame
    % For reference please refer the attached documents

    % Calculation of angular momentum
    ang_mom = cross(position,velocity);
    
    % Calculation of inclination
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.58
    if ang_mom(3)>0
        incli=atan(sqrt(ang_mom(1)^2+ang_mom(2)^2)/ang_mom(3));
    else
        incli=pi+atan(sqrt(ang_mom(1)^2+ang_mom(2)^2)/ang_mom(3));
    end
    incli=incli*180/pi;

    % Calculation of The Right ascension of the ascending node
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.58
    if ang_mom(2)<0
        if(ang_mom(1)>0)
            acend_node=atan(ang_mom(1)/-ang_mom(2));
        else
            acend_node=2*pi+atan(ang_mom(1)/-ang_mom(2));
            
        end
    else
        acend_node=pi+atan(ang_mom(1)/-ang_mom(2));
    end
    acend_node=acend_node*180/pi;

    % Calculation of semi latus rectum
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.59
    sem_latus = (norm(ang_mom))^2/(Grav_C*M_earth);

    % Calculation of Semi major axis
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.60
    semi_axis = 1/(2/norm(position)-(norm(velocity)^2)/(Grav_C*M_earth));

    % Calculatio of Mean motion
    % For reference refer book Satellite Orbits Page 28 ,Equation 2.61
    mean_m = sqrt(Grav_C*M_earth/semi_axis^3); 

    % Calculation of Eccentricity
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.62
    eccen = sqrt(1-(sem_latus/semi_axis));

    % Calculation of Eccentric Anomaly
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.64
    if (1-norm(position)/semi_axis)>0
        ecc_anom = atan(dot(position,velocity)/(semi_axis^2*mean_m*(1-norm(position)/semi_axis)));
    else
        ecc_anom = pi+ atan(dot(position,velocity)/(semi_axis^2*mean_m*(1-norm(position)/semi_axis)));
    end

    % Calculation of Mean Anomaly
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.65
	mean_anom = ecc_anom-eccen*sin(ecc_anom);
    if mean_anom>=0
        mean_anom=mean_anom*180/pi;
    else
        mean_anom=mean_anom*180/pi+360;        
    end
    
    % Calculation of Argument of perigee
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.66
    % Calculation of Unit vector Along Anglar momemtum Vector
    ang_mom_u=ang_mom/norm(ang_mom);
    % Calculation of Argument of Latitude
    if (-position(1)*ang_mom_u(2)+position(2)*ang_mom_u(1))>0
        arg_lat = atan(position(3)/(-position(1)*ang_mom_u(2)+position(2)*ang_mom_u(1)));
    else
        arg_lat =pi+ atan(position(3)/(-position(1)*ang_mom_u(2)+position(2)*ang_mom_u(1)));
    end

    % Calculation of True Anomaly
    % For reference refer book Satellite Orbits Page 29 ,Equation 2.67
    if (cos(ecc_anom)-eccen)>0
        true_anom = atan(sqrt(1-eccen^2)*sin(ecc_anom)/(cos(ecc_anom)-eccen));
    else
        true_anom=pi+atan(sqrt(1-eccen^2)*sin(ecc_anom)/(cos(ecc_anom)-eccen));
    end

    arg_per = arg_lat-true_anom;
    if (arg_per>0)
        arg_per = arg_per*180/pi;  
    else
        arg_per = arg_per*180/pi+360;
    end
    
    
end
    
    

    
    

