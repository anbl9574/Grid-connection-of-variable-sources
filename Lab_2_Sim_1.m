clear all
close all
clc
%% Loading Data (Do not Change)
folder = cd;
addpath([folder '\m-files'],[folder '\Radiation data']);
load('Barcelona.mat'); load('Stockholm.mat');

%% SIMULATION 1 - Yearly Irradiance

%%% Irradiation
Tilt = 53;           % Tilt of Solar Panel
Azimuth = 28;        % Azimuth of Solar Panel
Site = STOCKHOLM;    % Site of Solat Panel
Temperature = 05;    % Ambient Temperature

[IbT, IdT, IgT] = solrad(Site, Tilt, Azimuth, 0.2, 1);
Gt = IbT+ IdT+IgT;    % Calculate Global Irradiance throughout the year
Gt = Gt(1839)      % Calculate Global Irradiance In the specified hour

%%% Datasheet Solar Panel 250W
Tc_noct = 44; % Nominal Operation Cell Temperature (NOCT)
Ta_noct = 20; % Ambient Temperature for NOCT measurement
Gt_noct = 800; % Global Irradiance for NOCT measurement
n_ref   = 15.3/100; % Efficiency

%%% Cell Temperature Calculation
Ta = Temperature;
Tc = Ta+Gt*((Tc_noct-Ta_noct)/Gt_noct)*(1-n_ref)








