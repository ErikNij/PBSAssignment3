close all; 
clear;
clc
%% Initializing and getting correct file
current_directory = "C:\Users\20183809\Documents\GitHub\PBSAssignment3";

allvelo = fullfile(current_directory,"Transientvelocities.txt");
kbT = 1;
kb = 1.380649e-23;
m = 1;
binsize = 0.1;
xplot = -10:0.1:10;
data = load(allvelo);

%% Getting velocity and Maxwell boltzmann distribution
velocity = data((end-fix(height(data)*1/2)):end,:);
v_squared = mean(data.^2,"all");
Treal = m*mean(v_squared,"all");
Tboltz = m/kb*mean(velocity.*velocity,"all");
 
f = @(v)(sqrt(m/(2*pi*kb*Tboltz))*exp(-m*v.^2./(2*kb*Tboltz)));
[pfinal,edgesfinal] = histcounts(velocity,xplot);
edgesfinal(end) =[];

%% Figures
figure
subplot(3,1,1)
plot(edgesfinal,pfinal/trapz(binsize,pfinal));
title("Velocities of particles")
subplot(3,1,2)
plot(xplot,f(xplot)/trapz(binsize,f(xplot)));
string = sprintf("Temperature T in simulation: %f KbT units",Treal);
title(string)
subplot(3,1,3)
hold on
plot(edgesfinal,pfinal/trapz(binsize,pfinal));
plot(xplot,f(xplot)/trapz(binsize,f(xplot)));
hold off
