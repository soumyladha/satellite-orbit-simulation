%*************************************************************************
% This function creates a user interface, and calls other
% functions in order to calculate satellite pos, vel, acc and jerk
%  in the ECEF frame. The inputs are orbital parameters in ECEF  
%  frame at the input time or in J2000 ECI frame.
% Functions Called :  
%                   KepCartTime
%                   CartKep
% Library Calls :
%                   NIL
% Function Output :
%                   semi_axis = Semi major axis
%                   eccen = Eccentricty of Orbit
%                   incli = Inclination
%                   acend_node = The Right ascension of the ascending node 
%                                in degrees
%                   arg_per = Argument of perigee in degrees
%                   mean_anom = Mean anomaly in degrees
%                   pos_par = position
%                   vel_par = velocity
%                   acc_par = accelaration
%                   jerk_par = jerk
% Global Variables : NIL
% 
% 
% Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************
clc; % Clears command history
clearvars; % Clears workspace variables

% ===========================================================================
% In the message below, the flag for ECEF or ECI frame orbital parameters is
% missing. Also, UTC time of the first time step, duration of simulation and 
% time  step should be user inputs.
% When ECI is chosen, additional inputs are needed (UT1 - UTC in s, TT - UTC
% in s, x_p and y_p  for polar motion with their unit indicated)
% Both for precession and nutation the time is in TT. First, convert UTC to 
% Julian day and then add corrections for TT.
% ===========================================================================

message = ['The input parameters should be given in prescribed format::\n\n'...
    'When input is given from the file it should be separated \n'...
    'by space or semicolon\n\n'...
    'Semi major axis in (mtr)\n'...
    'Eccentricity of the orbit\n'...
    'Inclination (degrees)\n'...
    'Right ascension of the ascending node (degrees) \n'...
    'Argument of perigee (degrees)\n'...
    'Mean anomaly (degrees)\n'...
    'ECEF=0 / ECI =1\n'...
    'Date : format day/month/year/hours/minutes/seconds\n'...
    'Delta UT1:\n'...
    'Polar Motion x:\n'...
    'Polar Motion y:\n'...
    'TT-UTC \n'...
    'Step Size \n'...
    'Total Time \n'...
    '0 for manual INPUT\n1 for INPUT from FILE\n2 for default INPUT\n\n'];
fprintf(message) % Printing the message   
binary = input('Waiting for input:'); % Taking the input

if binary == 1
    
    % Exctracting user input from file
    % Extracting All except Julian date from input
    % Locating the Path and File Name
    [file_name,PathName] = uigetfile('*.txt','Select the INPUT file');
    format=('%s %f'); % Specifying format (string, float)
    fid=fopen(file_name,'r'); % File identifier
    cell_temp = textscan(fid,format); % Converting input to Cell
    input_orb_par = cell2mat(cell_temp(2)); % Converting input to matrix
    
   
elseif binary == 0
    
    % Extracting input from user input in text field
    prompt = {'Semi major axis (m):','Eccentricty of Orbit:','Inclination (degrees):',...
        'Right ascension of the ascending node (degrees):',...
        'Argument of perigee (degrees):','Mean/True anomaly (degrees):'...
        'flag = 0 ECEF & flag = 1 ECI:','Julian_date :',...
        'Delta UT1(degrees)','x:(degrees)','y:(degrees)','TT-UTC(degrees)','Step Size:','Total Time (hr):','Total Time (min):','Total Time (sec):',};
    dlg_title = 'Input'; % Heading of the dialog box
    num_lines = 1; % No of lines in the text box
    % Setting default values
    defaultans = {'2.655970469016551e+07','0.1836776733E-002','55.1402','-1.198805165047906e+02 ','-72.015595427841006','-53.573906420897934','1','1999,03,04,0,00,-13','0.649232','0.06740','0.24173','64.184','100','24','0','0'}; 
    input_orb_par = inputdlg(prompt,dlg_title,[1 50; 1 50; 1 50;1 50;1 50;1 50;1 50;1 50;...
        1 50;1 50;1 50;1 50;1 50;1 50;1 50;1 50;]...
        ,defaultans);
    temp_var = juliandate(str2num(cell2mat(input_orb_par(8))));
    input_orb_par = str2double(input_orb_par);
    input_orb_par(8) = temp_var; 
    
    button = questdlg('Do you want save as default?',...
    'Continue Operation','Yes','No','No');

    if strcmp(button,'Yes')
        % Writing the input values to the file
        disp('Creating file')
        fileID = fopen('newdefault.txt','w');
        fprintf(fileID,'%6s %12.8f\r\n','Semi_major_axis_(m)',input_orb_par(1));
        fprintf(fileID,'%6s %12.8f\r\n','Eccentricty_of_Orbit',input_orb_par(2));
        fprintf(fileID,'%6s %12.8f\r\n','Inclination_(degrees)',input_orb_par(3));
        fprintf(fileID,'%6s %12.8f\r\n','Right_ascension_of_the_ascending_node_(deg)'...
            ,input_orb_par(4));
        fprintf(fileID,'%6s %12.8f\r\n','Argument_of_perigee_(degrees):',input_orb_par(5));
        fprintf(fileID,'%6s %12.8f\r\n','Mean/True_anomaly_(degrees):',input_orb_par(6));
        fprintf(fileID,'%6s %12.8f\r\n','flag=0_ECEF_&_flag=1_ECI:',input_orb_par(7));
        fprintf(fileID,'%6s %12.8f\r\n','Julian_date',input_orb_par(8));
        fprintf(fileID,'%6s %12.8f\r\n','Delta_UT1',input_orb_par(9));
        fprintf(fileID,'%6s %12.8f\r\n','x_polar',input_orb_par(10));
        fprintf(fileID,'%6s %12.8f\r\n','y_polar',input_orb_par(11));
        fprintf(fileID,'%6s %12.8f\r\n','TT-UTC',input_orb_par(12));
        fprintf(fileID,'%6s %12.8f\r\n','Step_size',input_orb_par(13));
        fprintf(fileID,'%6s %12.8f\r\n','total_time(hr)',input_orb_par(14));
        fprintf(fileID,'%6s %12.8f\r\n','total_time(min)',input_orb_par(15));
        fprintf(fileID,'%6s %12.8f\r\n','total_time(sec)',input_orb_par(16));
        
        
    end
    

elseif binary == 2
    
    % Extracting the user input from new defaults file
    file=('newdefault.txt');
    format=('%s %f');
    fid=fopen(file);
    cell_temp = textscan(fid,format);
    input_orb_par = cell2mat(cell_temp(2));
end

%Degree to radian conversion of all angular parameters
for index = 3:6
    input_orb_par(index) = input_orb_par(index)*pi/180;
end

disp('Simulating')

[pos_par,vel_par,acc_par,jerk_par] = KepCartTime(input_orb_par(1),input_orb_par(2),...
    input_orb_par(3),input_orb_par(4),input_orb_par(5),input_orb_par(6),input_orb_par(7),...
    input_orb_par(8),input_orb_par(9),input_orb_par(10),input_orb_par(11),....
    input_orb_par(12),input_orb_par(13),3600*input_orb_par(14)+60*input_orb_par(15)+input_orb_par(16));

% Back calculation of orbital parameters

if input_orb_par(7)==0
    disp('Back Calculating')

    BackCalPar=zeros(size(pos_par,1),6);

    for index = 1:1:size(pos_par,1)
        [BackCalPar(index,1), BackCalPar(index,2), BackCalPar(index,3), BackCalPar(index,4), ...
        BackCalPar(index,5), BackCalPar(index,6)]=...
                            CartKep(pos_par(index,:),vel_par(index,:),input_orb_par(7),...
        input_orb_par(8));
    end
end
% For generating the dialog box
% Option for saving the result to file
button = questdlg('Do you want save to file',...
    'Continue Operation','Yes','No','No');
    if strcmp(button,'Yes')
        disp('Creating file')
        fileID = fopen('position.txt','w');
        fprintf(fileID,'%s %s %s\r\n','position(x)','position(y)','position(z)');
        fprintf(fileID,'%f %f %f\r\n',pos_par');
        fileID = fopen('velocity.txt','w');
        fprintf(fileID,'%s %s %s\r\n','velocity(x)','velocity(y)','velocity(z)');
        fprintf(fileID,'%f %f %f\r\n',vel_par');
        fileID = fopen('accelaration.txt','w');
        fprintf(fileID,'%s %s %s\r\n','accelaration','accelaration','accelaration');
        fprintf(fileID,'%f %f %f\r\n',acc_par');
        fileID = fopen('jerk.txt','w');
        fprintf(fileID,'%s %s %s\r\n','jerk','jerk','jerk');
        fprintf(fileID,'%f %f %f\r\n',jerk_par');
    end