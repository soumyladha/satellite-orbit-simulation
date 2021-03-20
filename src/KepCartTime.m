%*************************************************************************
% This function calls the function KepCart in order to calculate
% the position, velocity, accelaration and jerk of the satellite at different
% time by varying mean anomaly
%
% Functions called :
%                   KepCart
% Return Values :
%                   pos_vec = vector of positions 
%                   vel_vec = vector of velocity 
%                   acc_vec = vector of accelaration 
%                   jerk_vec = vector of jerk 
%
% Function Arguments :
%                   semi_axis = Semi major axis
%                   eccen = Eccentricty of Orbit
%                   incli = Inclination
%                   acend_node = The Right ascension of the ascending node in degree
%                   arg_per = argument of perigee in degree
%                   mean_anom = Mean anomaly in degree
%                   flag = it determines whether the output is in ECEF frame
%                   flag = 1 for ECI ;0 for ECEF
%                   julian_date = value of date in Julian
%  
% Library Calls :                  
%                   sqrt (square root) 
%                   Plot3 (For plotting 3-Dimensional graphs)
%                   disp (Shows display messages)
%                   prompt (user input Library)
%                   geoshow (show geographical map of earth)
%
%Function Output :  pos_vec
%		    vel_vec
%                   acc_vec
%                   jerk_vec
%
% Global Variables : NIL
% 
% 
% Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************

function [pos_vec,vel_vec,acc_vec,jer_vec]= KepCartTime(semi_axis,eccen,...
    incli,acend_node,arg_per,mean_anom,flag,julian_date,delta_ut1,x_polar,y_polar...
    ,tt_utc,step_size,total)

if flag == 1 % The orbital elements are in J2000 ECI frame
    %Loading Constants
    Constants;


    % Calculation of Mean Motion
    % For reference refer to book Satellite Orbits Page 23 ,Equation 2.35
    mean_m=sqrt(Grav_C*M_earth/semi_axis^3);

    % Calculation of time period

%      time_p = 2*pi/mean_m;
    index = 0; % index initialization
    
    
    
    % Preallocation
    vel_vec = zeros(ceil(total/step_size),3);
    pos_vec = zeros(ceil(total/step_size),3);
    acc_vec = zeros(ceil(total/step_size),3);
    jer_vec = zeros(ceil(total/step_size),3);
    latitude = zeros(ceil(total/step_size),1);
    longitude = zeros(ceil(total/step_size),1);
    
    for time = 0:step_size:total
        index = index + 1;
        
        % mean anomaly increament for each time 'T'
        mean_anom_temp = mean_anom + mean_m*time;
        [position_vec,velocity_vec,accion_vec,jerk_vec] = KepCart(semi_axis,eccen,incli,...
                    acend_node,arg_per,mean_anom_temp,flag,julian_date);
        
        %Conversion from ECI to ECEF at current time  
            
            
            coor_mat_eci_to_ecef_dt0 = ECI2ECEF(julian_date,delta_ut1,x_polar,y_polar,tt_utc);
            
            coor_mat_eci_to_ecef_dt1 = ECI2ECEF(julian_date+1/86400,delta_ut1,x_polar,y_polar,tt_utc);
            
            coor_mat_eci_to_ecef_dt2 = ECI2ECEF(julian_date+2/86400,delta_ut1,x_polar,y_polar,tt_utc);
            
            coor_mat_eci_to_ecef_dt3 = ECI2ECEF(julian_date+3/86400,delta_ut1,x_polar,y_polar,tt_utc);
            
            
            
            position_vec = position_vec*coor_mat_eci_to_ecef_dt0';

            
            % For Reference frame derivative--> 1st derivative
            ecef_2_eci_d1 = (coor_mat_eci_to_ecef_dt1'-coor_mat_eci_to_ecef_dt0')';
            velocity_vec = velocity_vec*coor_mat_eci_to_ecef_dt0'-position_vec*ecef_2_eci_d1*coor_mat_eci_to_ecef_dt0';
            
            
            % For reference frame derivatve-->  2nd derivative
            ecef_2_eci_d2 = (-2*coor_mat_eci_to_ecef_dt1'+...
                coor_mat_eci_to_ecef_dt0'+coor_mat_eci_to_ecef_dt2')';
            accion_vec = (accion_vec-...
                2*velocity_vec*ecef_2_eci_d1-position_vec*ecef_2_eci_d2)*coor_mat_eci_to_ecef_dt0';
            
            % For reference frame derivative-->  3rd derivative
            ecef_2_eci_d3 = (coor_mat_eci_to_ecef_dt3'-3*coor_mat_eci_to_ecef_dt2'...
                +3*coor_mat_eci_to_ecef_dt1'-coor_mat_eci_to_ecef_dt0')'*1;
            jerk_vec = (jerk_vec-3*accion_vec*ecef_2_eci_d1-3*velocity_vec*ecef_2_eci_d2-...
                position_vec*ecef_2_eci_d3)*coor_mat_eci_to_ecef_dt0';
            
  
        [latitude(index),longitude(index),~]=Geodetic(position_vec(1),position_vec(2),position_vec(3));
        % Save the outputs in arrays
        pos_vec(index,:) = position_vec;
        vel_vec(index,:) = velocity_vec;
        acc_vec(index,:) = accion_vec;
        jer_vec(index,:) = jerk_vec;
        
        julian_date=julian_date+step_size/(86400);
    
    end
        newpos=sqrt(sum(pos_vec.^2, 2));
        newvel=sqrt(sum(vel_vec.^2, 2));
        newacc=sqrt(sum(acc_vec.^2, 2));
        newjerk=sqrt(sum(jerk_vec.^2, 2));

figure(1)
        time = 0:step_size/3600:total/3600;
        subplot(2,2,3), hold on; plot(time,acc_vec(:,1),'r')
        subplot(2,2,3),hold on; plot(time,acc_vec(:,2),'g')
        subplot(2,2,3),hold on; plot(time,acc_vec(:,3),'b')
        subplot(2,2,3),hold on; plot(time,newacc(:,1),'k')        
        legend('x acceleration','y acceleration','z acceleration','Magnitude')
        ylabel('acceleration (mtr/sec^2)')
        xlabel('time (hr)')
        grid on;
        
figure(1)
        time = 0:step_size/3600:total/3600;
        subplot(2,2,1), hold on; plot(time,pos_vec(:,1),'r')
        subplot(2,2,1),hold on; plot(time,pos_vec(:,2),'g')
        subplot(2,2,1),hold on; plot(time,pos_vec(:,3),'b')
        subplot(2,2,1),hold on; plot(time,newpos(:,1),'k')
        legend('x position','y position','z position','Magnitude')
        ylabel('position (mtr)')
        xlabel('time (hr)')
        grid on;





figure(1)
        time = 0:step_size/3600:total/3600;
        subplot(2,2,2), hold on; plot(time,vel_vec(:,1),'r')
        subplot(2,2,2),hold on; plot(time,vel_vec(:,2),'g')
        subplot(2,2,2),hold on; plot(time,vel_vec(:,3),'b')
        subplot(2,2,2),hold on; plot(time,newvel(:,1),'k')
        legend('x velocity','y velocity','z velocity','Magnitude')
        ylabel('velocity (mtr/sec)')
        xlabel('time (hr)')
        grid on;


figure(1)
        time = 0:step_size/3600:total/3600;
        subplot(2,2,4), hold on; plot(time,jerk_vec(:,1),'r')
        subplot(2,2,4),hold on; plot(time,jerk_vec(:,2),'g')
        subplot(2,2,4),hold on; plot(time,jerk_vec(:,3),'b')
        subplot(2,2,4),hold on; plot(time,newjerk(:,1),'k')        
        legend('x jerk','y jerk','z jerk','Magnitude')
        ylabel('jerk (mts/sec^3)')
        xlabel('time (hr)')
        grid on; 
          
        
else % The orbital elements are in ECEF at the input time
    %Loading Constants
    Constants;


    % Calculation of Mean Motion
    % For reference refer book Satellite Orbits Page 23 ,Equation 2.35
    mean_m=sqrt(Grav_C*M_earth/semi_axis^3);
     

    % Calculation of time period
    time_p = 2*pi/mean_m;



    index = 0; % index initialization
    
    % Preallocation
    latitude = zeros(ceil(total/step_size),1);
    longitude = zeros(ceil(total/step_size),1);
    vel_vec = zeros(ceil(total/step_size),3);
    pos_vec = zeros(ceil(total/step_size),3);
    acc_vec = zeros(ceil(total/step_size),3);
    jer_vec = zeros(ceil(total/step_size),3);
    
    for time = 0:step_size:total
        index = index + 1;
        if (rem(index,floor((time_p/(step_size)))) < 1e-4)
            
            disp(strcat('orbit_no = ',int2str(ceil(time/time_p))));
        end
        
        % mean anomaly increament for each time 'T'
        mean_anom_temp = mean_anom + mean_m*time; % Without resetting
        
        % The Right ascension of the ascending node change for each time 'T'

        acend_node_temp = acend_node - norm(Omega_earth)*time; % Without resetting
        
        % Calculate pos, vel, acc and jerk for each time instant
        [position_vec,velocity_vec,accion_vec,jerk_vec] = KepCart(semi_axis,eccen,incli,...
                  acend_node_temp,arg_per,mean_anom_temp,flag);

    
        % Calculate the latitude and longitude and height from ECEF position
        
        
        % For reference please the refer the document...
        % GNSS,inertial and multisensor navigation by Paul D Groves 
        % Page no 38
        [latitude(index),longitude(index),~]=Geodetic(position_vec(1),position_vec(2),position_vec(3));

    
        % Save the outputs in arrays
        pos_vec(index,:) = position_vec;
        vel_vec(index,:) = velocity_vec;
        acc_vec(index,:) = accion_vec;
        jer_vec(index,:) = jerk_vec;
        
        
        
    
    end
        
    % Writing to the file
        fileID = fopen('latlog.txt','w');
        fprintf(fileID,'%s %s\r\n','Latitude','Longitude');
        fprintf(fileID,'%f %f\r\n',[latitude longitude]');
        
        newpos=sqrt(sum(pos_vec.^2, 2));
        newvel=sqrt(sum(vel_vec.^2, 2));
        newacc=sqrt(sum(acc_vec.^2, 2));
        newjerk=sqrt(sum(jer_vec.^2, 2));
       
        
figure(1)
        time = 0:step_size/3600:total/3600;
        subplot(2,2,1), hold on; plot(time,pos_vec(:,1),'r')
        subplot(2,2,1),hold on; plot(time,pos_vec(:,2),'g')
        subplot(2,2,1),hold on; plot(time,pos_vec(:,3),'b')
        subplot(2,2,1),hold on; plot(time,newpos(:,1),'y')
        legend('x position','y position','z position')
        ylabel('position (mtr)')
        xlabel('time (sec)')
        grid on;

figure(1)
        time = 0:step_size/3600:total/3600;
        subplot(2,2,2), hold on; plot(time,vel_vec(:,1),'r')
        subplot(2,2,2),hold on; plot(time,vel_vec(:,2),'g')
        subplot(2,2,2),hold on; plot(time,vel_vec(:,3),'b')
        subplot(2,2,2),hold on; plot(time,newvel(:,1),'y')
        legend('x velocity','y velocity','z velocity')
        ylabel('velocity (mtr/sec)')
        xlabel('time (hr)')
        grid on;
        
figure(1)
        time = 0:step_size/3600:total/3600;
        subplot(2,2,3), hold on; plot(time,acc_vec(:,1),'r')
        subplot(2,2,3),hold on; plot(time,acc_vec(:,2),'g')
        subplot(2,2,3),hold on; plot(time,acc_vec(:,3),'b')
        subplot(2,2,3),hold on; plot(time,newacc(:,1),'y')        
        legend('x acceleration','y acceleration','z acceleration')
        ylabel('acceleration (mtr/sec^2)')
        xlabel('time')
        grid on;
        
figure(1)
        time = 0:step_size/3600:total/3600;
        subplot(2,2,4), hold on; plot(time,jer_vec(:,1),'r')
        subplot(2,2,4),hold on; plot(time,jer_vec(:,2),'g')
        subplot(2,2,4),hold on; plot(time,jer_vec(:,3),'b')
        subplot(2,2,4),hold on; plot(time,newjerk(:,1),'y')        
        legend('x jerk','y jerk','z jerk')
        ylabel('jerk (mtr/sec^3)')
        xlabel('time (hr)')
        grid on;
        
figure(2)
        %Plotting ground track of the satellite
        
        disp('plotting Ground Track')
        grid on
        geoshow('landareas.shp', 'FaceColor', [1 1 1])% Plotting world map
        
        axis off
        axis([-200 200 -90 90]);% Limiting the scale of axis
        axis on
        hold on
        grid on;
        box on
        plot(longitude,latitude,'b');
        hold off    
        
      
        
end

