%*************************************************************************
% This function does the rotation about z axis with counter clockwise as
% positive direction 
% Function Argument :
%                     z (angle in radian)
% Library calls : NIL
% Functions calls : NIL
% Global Variables : NIL
% Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************

function rotz = Rotz(z)
% For reference refer to GNSS,inertial and multisensor navigation...
% by Paul D Groves, Page 26

rotz=[cos(z) sin(z) 0;-sin(z) cos(z) 0;0 0 1];
end