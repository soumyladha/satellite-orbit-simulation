%*************************************************************************
% This function does the rotation about y axis with counter clockwise as
% positive direction 
% Function Argument :
%                     y (angle in radian)
% Library calls : NIL
% Functions calls : NIL
% Global Variables : NIL
% Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************

function roty = Roty(y)
% For reference refer to GNSS,inertial and multisensor navigation...
% by Paul D Groves, Page 26

roty=[cos(y) 0 -sin(y);0 1 0;sin(y) 0 cos(y)];
end