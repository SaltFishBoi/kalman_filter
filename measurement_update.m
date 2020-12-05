% measurement_update.m
function [x,P]=measurement_update(z, H, xm, Pm, R) 

        K=Pm*(H')/(H*Pm*(H')+R); 
        %New state 
        x=xm+K*(z-H*xm); 
        P=(eye(size(Pm))-K*H)*Pm; 

end 