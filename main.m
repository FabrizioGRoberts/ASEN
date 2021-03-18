%Fabrizio Roberts
%Student ID: 109582514
%Last Edit: December 4th, 2020
%2012 Project 2

%This code etimates the trajectory of a bottle rocket given the parameters
%below with the use of ODE45

clc
clear
close all

%constant variables (2D: [x z])
global CD Cdis At Ab dT dB g gam pA pW Vem R stL Pat Vair0 P0 mair0
tspan = [0 5]; %time span (s)
stL = 0.5; %length of test stand (m)
CD = 0.5; %drag coefficient
Cdis = 0.8; %discharge coefficient
Vem = 0.002; %volume of empty bottle (m^3)
Pat = 83426.56; %atmospheric pressure [12.1 psi](pa)
Ti = 300; %initial air temp (K)
gam = 1.4; %ratio of specific heats for air
dT = 0.021; %diameter of bottle throat (m)
At = pi*((dT/2)^2); %cross sectional area of bottle throat (m^2)
dB = 0.105; %diameter of bottle 
Ab = pi*((dB/2)^2); %cross sectional area of bottle (m^2)
g = -9.81; %gravity vector (m/s^2)
R = 287; %air gas constant (J/kgK)
P0 = 428164.56; %initial pressure in bottle(pa)
mB = 0.15; %mass of empty bottle (kg)
pA = 0.961; %air density (kg/m^3)
pW = 1000; %water density (kg/m^3)
Vw = 0.001; %initial volume of water (m^3)
Vair0 = Vem-Vw; %initial air volume

%initial values
v0 = 0; %initial velocity (m/s)
x0 = 0; %initial x position (m)
z0 = 0.25; %initial z position (m)
theta0 = pi/4; %initial angle of rocket (radians)
mair0 = (P0*Vair0)/(R*Ti); %air mass initial (kg)
mR0 = mair0+mB+(Vw*pW); %rocket mass initial (kg)


%ODE45
X0 = [x0; z0; v0; v0; theta0; mR0; mair0; Vair0];
st = odeset('Events',@stop);
[t,X] = ode45(@(t,X) rocket(t,X),tspan,X0,st);

%Plotting Thrust
figure
[TF,Phase] = TF(X);
plot(t,TF);
title("Thrust")
xline(t(43),'-.',"Phase 2");
xline(t(66),'-.',"Phase 3");
xlim([0 0.5])
ylim([0 300])
xlabel("Time (s)");
ylabel("Thrust (N)");

%Plotting Trajecttory
figure
plot(X(:,1),X(:,2));
xlabel("X Position (m)");
ylabel("Z Position (m)");
title("Bottle Rocket Trajectory")
xl = xline(X(43,1),'-.',"Phase 2");
xl.LabelHorizontalAlignment = 'left';
xline(X(66,1),'-.',"Phase 3");
xlim([0 100]);
ylim([0 30]);

%Plotting Z Component of Velocity
figure
plot(t,X(:,4))
title("Z Velocity")
xlim([0 4]);
ylim([-20 25]);
xlabel("Time (s)");
ylabel("Z-Velocity (m/s)");

%Plotting X Component of Velocity
figure
plot(t,X(:,3))
title("X Velocity")
xlim([0 4]);
ylim([0 30]);
xlabel("Time (s)");
ylabel("X-Velocity (m/s)");

%Plotting Volume
figure
plot(t,X(:,8))
title("Volume")
xlim([0 0.25])
ylim([(1*10^-3) (2*10^-3)])
xlabel("Time (s)");
ylabel("Volume of Air in Bottle Rocket (m^3)");
