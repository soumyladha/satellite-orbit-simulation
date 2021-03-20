%*************************************************************************
% This function does the rotation about x axis with counter clockwise as
% positive direction 
% Function Argument :
%                     x (angle in radian)
% Library calls : NIL
% Functions calls : NIL
% Global Variables : NIL
% Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************

function rotx = Rotx(x)
% For reference refer to GNSS,inertial and multisensor navigation...
% by Paul D Groves, Page 26
rotx=[1 0 0;0 cos(x) sin(x);0 -sin(x) cos(x)];
end