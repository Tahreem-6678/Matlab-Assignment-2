clear all
%% Constants
R_E=6371.2e3;                                                 % Earth's radius [m]
mu=3.986e14;                                                  % Earth's gravitational parameter [m^3*s^-2]
mu_0=1.25663706212e-6;                                        % Vacuum permeability [H/m]
Earth_period=86164;                                           % Siderial period of Earth [s]
Earth_ang_speed=2*pi/Earth_period;                            % Earth's angular speed

%% Satellite specifications (mockup)
%dim_1=0.1;                                                   % Dimension 1 [m]
%dim_2=0.1;                                                    % Dimension 2 [m]
%dim_3=0.3;                                                    % Dimension 3 [m]
%sat_mass=150;                                                   % Mass[Kg]
%J_11=sat_mass/12*(dim_2^2+dim_3^2)-0.005;                           % Moment of inertia along 1 [Kg*m^2]
%J_22=sat_mass/12*(dim_1^2+dim_3^2)-0.001;                           % Moment of inertia along 2 [Kg*m^2]
%J_33=sat_mass/12*(dim_1^2+dim_2^2)+0.006;                           % Moment of inertia along 3 [Kg*m^2]

%matrix of inertia of my mockup already calculated
J_11=0.0077; J_22=0.0081; J_33=0.0069;                      %  𝐼xx =19.5  kg.m^2 ,𝐼𝑦𝑦=19 kg.m^2 ,𝐼𝑧𝑧=12.6 kg.m^2
J = [J_11  0   0                                              % Matrix of inertia
    0  J_22   0
    0  0   J_33];
invJ=inv(J);

%% Reaction wheel specifications
wheel_speed_max = 8000*2*pi/60;                                % Wheel maximum speed       [rad/s]
wheel_torque_max= 0.05;                                    % Maximum wheel torque      [N*m]
wheel_momentum_max=1.77e-3;                                   % Wheel maximum momentum storred [N*m*s]
J_w=wheel_momentum_max/wheel_speed_max;                       % Wheel moment of inertia     [Kg*m^2]
K_T=3.5e-3;                                                   % Torque constant     [N*m/A]
J_w_set=J_w*eye(3);
invJ_J_w=inv(J-J_w_set);

%% Magnetorquers specification
max_dipole=20;                                                 % Maximum dipole A*m^2

%% Orbit
ecc=0;                                                        % Eccentricity e
RAAN=0*pi/180;                                                % Rigth ascension of ascending node Omega [rad]
orb_i=0*pi/180;                                              % Inclination will need updating [rad]
arg_per=0*pi/180;                                             % Argument of perigee [rad]
altitude=400e3;                                               % Altitude [m]
sma=R_E+altitude;                                             % Semi-major axis a [m]
orb_ang_vel=sqrt(mu/sma^3);                                   % Orbital angular velocity or mean motion n [rad/s]
orb_period=2*pi/orb_ang_vel;                                  % Orbital period [s]

%% INITIAL CONDITION
yaw_0=0*pi/180;                                            % Initial rotation around Z body axis [rad]
pitch_0=0*pi/180;                                              % Initial rotation around Y [rad]
roll_0=0*pi/180;                                               % Initial rotation around X [rad]
DCM_0=angle2dcm(yaw_0,pitch_0,roll_0);                        % Initial DCM
q_0=dcm2quat(DCM_0);                                          % Initial quaternion
% w_0=[-0.0063   -0.0026    0.0025];
w_0 = [0  0 0];

%% Geomagnetic field model - IGRF 13 (2020) %%
% % Centered dipole model
% g_1_0=-29404.8;                                               % Degree n=1 order m=0
% g_1_1=-1450.9;                                                % Degree n=1 order m=1
% h_1_1=4652.5;                                                 % Degree n=1 order m=1
% B_0=sqrt(g_1_0^2+g_1_1^2+h_1_1^2);                            % Magnetic field intensity at the Earth's surface [nT]
% c_1_1=sqrt(g_1_1^2+h_1_1^2);
% theta_n=atan2(c_1_1/B_0,-g_1_0/B_0);                          % Colatitude magnetic dipole (ECEF) [rad]
% phi_n=atan2(-h_1_1/c_1_1,-g_1_1/c_1_1);                       % Longitude magnetic dipole (ECEF) [rad]
% M=4*pi/mu_0*B_0*(sma/1000)^3;                                 % Magnetic dipole intensity [A*m^2]

%% Control 
k_q = 0.0001;
k_w = 0.001;