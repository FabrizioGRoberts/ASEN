%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                               %
%                                                               %
%                       ASEN 2003 -- Lab 1                      %
%                        Fabrizio Roberts                       %
%           Roller Coaster Design Analysis -- Main Code         %
%     Created: June 01, 2021 -- Last Modified: June 08, 2021    %
%                                                               %
%                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
close all
clear
clc

%% Constants + Limits + Initial Values
g = 9.81;            %Force of gravity [m/s^2]            
h0 = 125;           %Starting height [m]
m = 500;            %Loaded vehicle mass [kg]
track0 = 1250;      %Total track length [m]
pos0 = [0 0 125];   %Starting position [x y z]
G0 = [1 0 0];       %Starting G-forces [up/down front/back left/right]
vel0 = 0;           %Starting velocity [m/s]
Gmax_up = 6;        %Max upward Gs allowed
Gmax_lateral = 3;   %Max lateral Gs allowed
Gmax_back = 4;      %Max braking Gs allowed

%Initialize Coaster State
state = [pos0 vel0 G0 track0];  %[xyz vel G-forces track-remaining]

%% Coaster Path Function Calls

%Free-fall Portion
[state_fall,d1] = fall(pos0,vel0,g,h0);                                              %Calls fall function

%Setting end-state values from free-fall portion
pos_fall_end = [state_fall(end,1) state_fall(end,2) state_fall(end,3)];         %Position end-state
vel_fall_end = state_fall(end,4);                                               %Velocity end-state
G_fall_end = state_fall(end,5);                                                 %G-Force end-state
r_valley = 15;                                                                  %Set circle radius

%Circular Valley Portion
[state_valley,d2] = valley(pos_fall_end, vel_fall_end, r_valley, G_fall_end, g, h0); %Calls valley function

%Setting end-state values from circular valley portion
pos_valley_end = [state_valley(end,1) state_valley(end,2) state_valley(end,3)]; %Position end-state
vel_valley_end = state_valley(end,4);                                           %Velocity end-state
G_valley_end = state_valley(end,5);                                             %G-Force end-state

%Parabolic Hill Portion
[state_hill,d3] = hill(pos_valley_end, vel_valley_end, G_valley_end, g, h0);         %Calls hill function

%Setting end-state values from parabolic hill portion
pos_hill_end = [state_hill(end,1) state_hill(end,2) state_hill(end,3)];
vel_hill_end = state_hill(end,4);
G_hill_end = state_hill(end,5);
theta_hill_end = atan(pos_hill_end(3)/pos_hill_end(1));
vel_hill_end2 = [vel_hill_end*cos(theta_hill_end), 0, vel_hill_end*sin(theta_hill_end)];
%%Loop Portion
[xyz, vf, vel_loop, g_force,d4] = loop2(pos_hill_end, vel_hill_end2, 22, 1.1064);
for i = 1:100
   gFl(i) = norm(g_force(i,:));
end
state_loop = [xyz vel_loop' gFl'];

%Banked Turn
[xyzb, vfb, vel_bank, g_forceb,d5] = banked(xyz(end,:),vf,8,0.5);
for i = 1:100
   gF(i) = norm(g_forceb(i,:));
end
state_bank = [xyzb vel_bank' gF'];

% %Dip to z = 0
% theta = pi/4;
% [xyzd, vd, g_forced] = dip(xyzb(end,:),vel_bank(end),theta,gF(end),25);
% state_dip = [xyzd, vd', g_forced'];

[xyzf, vel_out,velfin, g_forcef,d6] = parabolicHillSec(xyzb(end,:), vfb', 3.25);
state_fin = [xyzf, velfin', g_forcef'];

[xyzS,v,g_forceS,d7] = braking(xyzf(end,:),velfin(end),100);

state_stop = [xyzS, v', g_forceS(:,3)];
xG = zeros(1,700);
frtG = [xG g_forceS(:,1)'];
% Assembling Data
dTrav = d1+d2+d3+d4+d5+d6+d7;
state_final = cat(1,state_fall,state_valley,state_hill,state_loop,state_bank,state_fin,state_stop);    %Concatenate the track element states
%state_final = cat(1,state_fall,state_valley,state_hill,state_loop,state_bank,state_dip,state_stop);
%% Plotting

%Full coaster path plot w/ colored velocity
figure;
fig0 = plot3(state_final(:,1),state_final(:,2),state_final(:,3));

title('Rollercoaster Full Path')
grid on
colormap('turbo')
state_final(end,4) = NaN;
C = state_final(:,4);
patch(state_final(:,1),state_final(:,2),state_final(:,3),C,'EdgeColor','interp','FaceColor','none');
c2 = colorbar;
c2.Label.String = 'Velocity [m/s]';
xlabel('x position [m]')
ylabel('y position [m]')
zlabel('z position [m]') 
figure;

%G-Force plots
gplot = tiledlayout(3,1);
title(gplot,'G-forces vs. Track Location')

xstate = (1:length(state_final(:,1)));

nexttile
plot(xstate',state_final(:,5))
title('Up/Down G-forces')
xlabel('Track Piece (of 800)')
ylabel('G-Force')

nexttile
plot(xstate',zeros(1,800))
title('Lateral G-forces')
xlabel('Track Piece (of 800)')
ylabel('G-Force')

nexttile
plot(xstate',frtG)
title('Front/Back G-forces')
xlabel('Track Piece (of 800)')
ylabel('G-Force')



%% Track Element Functions

%Free-fall
function [state,d] = fall(pos,vel,g,h0)
     len_fall = 30;                                 %Length of free-fall [m]
     range = linspace(0,len_fall);                  %Set range of values for use
     
     x_fall = zeros(1,length(range));               %Preallocate x positions
     
     z_fall = zeros(1,length(range));               %Preallocate z positions
     z_fall(1) = pos(3);                            %Set iniital z position
     
     vel_fall = zeros(1,length(range));             %Preallocate velocity
     vel_fall(1) = vel;                             %Set initiial velocity
     
     G_fall = zeros(1,length(range));               %Preallocate G-forces
     
     %Loop to calculate height + velocity
     for j = 2:length(range)                
        z_fall(j) = h0 - range(j);
        vel_fall(j) = sqrt(2*g* (h0 - z_fall(j)));
     end
     
     %Set positions
     pos_fall(:,1) = x_fall';
     pos_fall(:,2) = 0;
     pos_fall(:,3) = z_fall';
     
     d = 30;
     
     %Generate state output
     state = [pos_fall vel_fall' G_fall'];
end

%Circular Valley
function [state,d] = valley(pos,vel,r,G,g,h0)
    theta_valley = linspace(0,(3*pi)/4);        %Theta range for element [rad]
    
    z_valley = zeros(1,length(theta_valley));   %Preallocate z values
    z_valley(1) = pos(3);                       %Set initial z value
    x_valley = zeros(1,length(theta_valley));   %Preallocate x values
    x_valley(1) = pos(1);                       %Set initial x value
    
    vel_valley = zeros(1,length(theta_valley)); %Preallocate velocity values
    vel_valley(1) = vel;                        %Set initial velocity
    
    G_valley = zeros(1,length(theta_valley));   %Preallocate G-force Values
    G_valley(1) = G(1);                         %Set initial G-force
    
    %Loop to calculate x,z positions, velocity and G-force
    for k = 2:length(theta_valley)
        x_valley(k) = x_valley(1) + (r - r*cos(theta_valley(k)));
        z_valley(k) = z_valley(1) - (r*sin(theta_valley(k)));
        vel_valley(k) = sqrt(2*g*(h0 - z_valley(k)));
        G_valley(k) = (vel_valley(k)^2 / (g*r)) - sin(theta_valley(k));
    end
    
    %Set Positions
    pos_valley(:,1) = x_valley';
    pos_valley(:,2) = 0;
    pos_valley(:,3) = z_valley';
    
    d = r*((3*pi)/4);
    
    %Generate state output
    state = [pos_valley vel_valley' G_valley'];
end

%Zero-G Hill
function [state,d] = hill(pos,vel,G,g,h0)
    theta0 = pi/4;                                  %Initial launch theta [rad]
    vel0 = vel;                                     %Initial velocity [m/s]
    vel0_z = vel0*sin(theta0);                      %Initial velocity in z direction [m/s]
    vel0_x = vel0*cos(theta0);                      %Initial velocity in x direction [m/s]
    
    range_hill = ((vel0^2*sin(2*theta0))/g);          %Calculate total horizontal range [m]
    x_range = linspace(pos(1),pos(1)+range_hill);   %Set points to calculate from along range
    z_hill = zeros(1,length(x_range));              %Preallocate z values
    z_hill(1) = pos(3);                             %Set initial z value
    vel_hill = zeros(1,length(x_range));            %Preallocate velocity values
    vel_hill(1) = vel;                              %Set initial velocity
    
    G_hill = zeros(1,length(x_range));              %Preallocate G-force values
    G_hill(1) = 0;                                  %Set initial G-force
    
    %Loop to calculate x,z positions, velocity and g-forces
    for m = 2:length(x_range)
        x_delta = x_range(m)- pos(1);
        z_hill(m) = z_hill(1) + (vel0_z*(x_delta/vel0_x)) - (.5*g*((x_delta/vel0_x)^2));
        vel_hill(m) = sqrt(2*g*(h0 - z_hill(m)));
        G_hill(m) = 0;
    end
    
    %Set positions
    pos_hill(:,1) = x_range';
    pos_hill(:,2) = 0;
    pos_hill(:,3) = z_hill';
   
    a = range_hill;
    b = x_range(end)-x_range(1);
    
    d = 0.5*sqrt((b^2)+16*(a^2))+((b^2)/(8*a))*log((4*a+sqrt((b^2)+16*(a^2)))/b);
    
    %Generate state output
    state = [pos_hill vel_hill' G_hill'];
end

%Loop
function [xyz, vf,vel_loop, g_force,d] = loop2(xyz0,v0,r,per)
%calculates theta in and, from that, builds a loop
% r = radius
% per = percent of circle wished to be completed

%% declaring variables
n = 100; %array size
h0 = 125; %m
g = 9.81; %m/s^2

%% calculating entry angle
x0 = xyz0(1);
y0 = xyz0(2);
z0 = xyz0(3);
v0_x = v0(1);
v0_z = v0(3);
theta0 = atan(v0_z/v0_x) + pi;

%% calcuating position matrix
if(v0_x < 0)
    theta = linspace(theta0,(2*pi*per+theta0),n);
else
    theta = linspace(theta0,(-2*pi*per+theta0),n);
end

x = x0 + r*sin(theta) - r*sin(theta0);
y = y0*ones(1,n);
z = z0 + r*cos(theta) - r*cos(theta0);

xyz = [x' y' z'];
vel_loop = (2*g*(h0 - z)).^0.5;

%% calculating final velocity
thetaf = 2*pi*per + theta0;
deltah = h0 - z(end);
spd = sqrt(2*g*deltah);

vf = [spd*cos(thetaf),0,spd*sin(thetaf)];

%% G force calculation
v = sqrt(2*g*(h0-z));
Gs = (v.^2)/(r*g) + cos(theta);
g_force = [zeros(n,1) zeros(n,1) Gs'];

d = 2*pi*r;
end

%Bank
function [xyz, vf, vel_bank, g_force,d] = banked(xyz0,v0,r,per)
%% Code for Banked Turn section
% r = radius
% per = percent of circle wished to be completed

%% declaring variables
n = 100; %array size
h0 = 125; %m
g = 9.81; %m/s^2
x0 = xyz0(1);
y0 = xyz0(2);
z0 = xyz0(3);
v0_x = v0(1);

%% calculating position matrix

%no lateral force function
thetaB = acot((r*g)/(v0_x^2)); %bank angle (radians)
theta = linspace(0,2*pi*per,n);
x = x0 + r*sin(theta);
y = y0 + r*cos(theta) - r;
z = z0*ones(1,n);
xyz = [x',y',z'];

%% calculating final velocity
thetaf = theta(n);
deltah = h0 - z(end);
spd = sqrt(2*g*deltah);
vf = [spd*cos(thetaf),0,spd*sin(thetaf)];
vel_bank = (2*g*(h0 - z)).^0.5;

%% G force calculation
Gs = (1/cos(thetaB))*ones(1,n);
g_force = [zeros(n,1),zeros(n,1),Gs'];

d = r*2*pi*per;
end

%Brakes
function [xyz,v,g_force,d] = braking(xyz0,v0,dis)
%Code for brake section

%% declaring variables
n = 100; %array size
g = 9.81; %m/s^2
x0 = xyz0(1);
y0 = xyz0(2);
z0 = xyz0(3);
v0_x = v0;

%% calculating required acceleration and Gs
a = -v0_x^2/(2*dis);
g_force = [(a/g)*ones(n,1),zeros(n,1),ones(n,1)];

%% calculating position matrix
x = linspace(0,dis,n);
vx = sqrt(v0_x^2+2*a*x);
if(v0_x>0)
    x = x0 - x;
else  
    x = x0 + x;
end
vf = [vx',zeros(n,1),zeros(n,1)];
xyz = [x',y0*ones(n,1),z0*ones(n,1)];
v = vx;
d = dis;
end

%parabolic hill part 2
function [xyz, vel_out, velfin, g_force,d] = parabolicHillSec(xyz_start, vel_in, t)
    %build path based on projectile motion calculations.
    %as the particle follows this path, it will have no outside forces
    %acting on it except for gravity.
    % grab all of the starting positions from the vector
    xin = xyz_start(1);
    yin = xyz_start(2);
    zin = xyz_start(3);
    n = 100;
    g = 9.81;
    h0 = 125;

        
    vx = -vel_in(1);
    vz = vel_in(3);
    % define t vector
    t = linspace(0,t);
    % calculate x
    x = xin - vx*t;
    z = zin + vz*t-.5*g*t.^2;

    % Compile the g forces and xyz coordinates into the matrices to be outputted
    g_forces = zeros(100,1);   % G-force matrix [front/back, left/right, up/down]

    vel_out = [vx, 0, vz-9.81*t(end) ];
    velfin = (2*g*(h0 - z)).^0.5; 

    % leveling section
    thetaIn = tan(z(end)/x(end));
    theta = linspace(thetaIn,(pi)/2);
    r = 40;
   
    z_valley = zeros(1,length(theta));   %Preallocate z values
    z_valley(1) = z(end);                       %Set initial z value
    x_valley = zeros(1,length(theta));   %Preallocate x values
    x_valley(1) = x(end);                       %Set initial x value
    
    vel_valley = zeros(1,length(theta)); %Preallocate velocity values
    vel_valley(1) = velfin(end);                        %Set initial velocity
    
    g_force = zeros(1,length(theta));   %Preallocate G-force Values
    g_force(1) = 0; 
      %Loop to calculate x,z positions, velocity and G-force
    for k = 2:length(theta)
        x_valley(k) = x_valley(1) - (r + r*cos(theta(k)));
        z_valley(k) = z_valley(1) - (r*sin(theta(k)))+12.2;
        vel_valley(k) = sqrt(2*g*(h0 - z_valley(k)));
        g_force(k) = (vel_valley(k)^2 / (g*r)) - sin(theta(k));
    end
    x = cat(1,x',x_valley')';
    y = yin*ones(1,length(x));
    z = cat(1,z',z_valley')';
    xyz = [x' y' z'];                           %XYZ matrix
    velfin = cat(1,velfin',vel_valley')';
    g_force = cat(1,g_forces,g_force')';
    
    a = zin-z_valley(1);
    b = 2*(xin-x_valley(1));
    
    dpa = (0.5*sqrt((b^2)+16*(a^2))+((b^2)/(8*a))*log((4*a+sqrt((b^2)+16*(a^2)))/b))/2;
    dva = r*(thetaIn - pi/2);
    d = dpa+dva;
end


