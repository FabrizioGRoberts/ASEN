function [dXdt] = rocket(t,X)
%ODE45 call: Caluculates change in X at each phase
%Input: State and Time (single step)
%Output: Change in State

global CD Cdis At Ab g gam pA pW Vem R stL Pat Vair0 P0 mair0

x = X(1); %x position
z = X(2); %z position
vx = X(3); %x velocity
vz = X(4); %z velocity
theta = X(5); %angle
mR = X(6); %rocket mass
mAir = X(7); %air mass
Vair = X(8); %current volume of air in rocket

v = sqrt((vx^2)+(vz^2));

vx = v*cos(theta);
vz = v*sin(theta);

Dx = 0.5*pA*(vx^2)*CD*Ab; %Drag x (N)
Dz = 0.5*pA*(vz^2)*CD*Ab; %Drag z (N)

if ~(Vair<Vem)
    Pend = P0*(Vair0/Vem)^gam; %Pressure 
    P2 = Pend*(mAir/mair0)^gam; %Internal pressure of rocket (phase 2)
    Vair = Vem;
end
if(Vair<Vem)
    %Phase 1:Water Propulsion
    P = P0*(Vair0/Vair)^gam; %air pressure function
    Ve = sqrt((2*(P-Pat))/pW); %exhaust velocity
    
    %calculating forces
    TF = 2*Cdis*At*(P-Pat); %Thrust Force vector of rocket (N)
    
    %variable rates
    mDot = -Cdis*At*Ve*pW; %mass flow rate
    mDotAir = 0; %air mass flow rate
    dVdt = Cdis*At*Ve; %rate of change in volume
elseif(Vair==Vem && P2>Pat)
    %Phase 2: Air Propulsion
    p = mAir/Vem;
    Pcrit = P2*((2/(gam+1))^(gam/(gam-1))); %Critical pressure
    T = P2/(R*p); %Air temperature
    
    %determining if flow choked
    if(Pcrit>Pat)
        %Exit Mach == 1
        Te = T*(2/(gam+1));
        pe = Pcrit/(R*Te);
        Ve = sqrt(gam*R*Te);
        Pe = Pcrit;
    else
        %Exit Mach != 1
        Me = sqrt((((P2/Pat)^((gam-1)/gam))-1)/((gam-1)/2));
        Te = T/(1+(((gam-1)/2)*(Me^2)));
        Ve = Me*sqrt(gam*R*Te);
        pe = Pat/(R*Te);
        Pe = Pat;
    end
    
    %variable rates
    mDot = -Cdis*pe*At*Ve; %mass flow rate
    mDotAir = mDot; %air mass flow rate
    TF = (Cdis*pe*At*Ve*Ve)+((Pat-Pe)*At); %Thrust Force vector of rocket (N)
    dVdt = 0; %change in volume
elseif(P2<=Pat)
    %Phase 3: Ballistic Flight
    TF = 0;
    mDot = 0;
    mDotAir = 0;
    dVdt = 0;
end

%acceleration calculated by component
ax = (TF*cos(theta) - Dx)/mR;
az = (TF*sin(theta) - Dz+(mR*g))/mR;

if(v<2 && x<(0.5/cos(theta))) %z<(stL*sin(pi/4))
    thetaDot = 0;
else
    thetaDot = (g*cos(theta))/v;
end

%returning change
dXdt(1) = vx; %x
dXdt(2) = vz; %z
dXdt(3) = ax; %x velocity {Issue}
dXdt(4) = az; %z velocity
dXdt(5) = thetaDot; %angle
dXdt(6) = mDot; %mass
dXdt(7) = mDotAir; %air mass
dXdt(8) = dVdt; %volume

dXdt = dXdt';
end

