close all
clear
clc
data = readtable("exp_1.txt");


% Extract the columns of interest (columns 6, 9, and 11)
time = data.Var1;
column6 = data.Var3;%(:, 6);
column8 = data.Var4;%(:, 8);
column10 = data.Var5;%(:, 10);
n = 600;
latest_points = column8(end - n + 1:end, 1);
%latest = table2array(latest_points);
mean = mean(latest_points);
% Create a figure for the plots
figure;


% Plot the data from column 6
subplot(3, 1, 1);
plot(time, column6)
%stackedplot(column6);
xlabel('Time [fs]');
ylabel('E [yg A^2 fs^{-2}]');
title('Potential energy');


% Plot the data from column 9
subplot(3, 1, 2);
plot(time, column8);
xlabel('Time [fs]');
ylabel('E [yg A^2 fs^{-2}]');
title('Kinetic energy');


% Plot the data from column 11
subplot(3, 1, 3);
plot(time, column10);
xlabel('Time [fs]');
ylabel('E [yg A^2 fs^{-2}]');
title('Total energy');