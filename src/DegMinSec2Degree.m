%*************************************************************************
% This function converts degree minutes seconds to degrees 
% Function Argument :
%                     degree 
%                     minutes
%                     seconds
% Library calls : NIL
% Functions calls : NIL
% Global Variables : NIL
% % Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************

function degrees = DegMinSec2Degree(degree,minutes,seconds)

    degrees = degree + minutes/60 + seconds/3600;

end