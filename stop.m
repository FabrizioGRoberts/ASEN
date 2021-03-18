function [zpos,stopCon,d] = stop(t,X)
%Stops the ODE45 call and z=0
    zpos = X(2)<0;
    stopCon = 1;
    d = 0;
end

