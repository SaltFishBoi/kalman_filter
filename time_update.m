%time_update.m

function [xm,Pm]=time_update(x, phi, P, Q) 
        xm=phi*x; 
        Pm=phi*P*(phi')+Q; 
end 