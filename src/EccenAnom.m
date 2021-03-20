%*************************************************************************
% The following function calculates the Eccentric anomaly from Mean Anomaly
% Function Arguments : ecc = Eccentricity
%                      mean_anom = Mean Anomaly
% Functions Outputs : Ecc_anom = Eccentric Anomaly
% Library calls : NIL
% Global Variables : NIL
% Version History: 
%                   <1.1> <Soumy Ladha>
%*************************************************************************

function ecc_anom = EccenAnom(ecc,mean_anom)

ecc_anom=0;
temp=mean_anom;
for i = 1:7
    % The method has been extracted ...
    % Page no 24 ,Eqn no 2.42
     ecc_anom=temp;
     temp=ecc_anom-(ecc_anom-ecc*sin(ecc_anom)-mean_anom)/(1-ecc*cos(ecc_anom));
end

