% main2D.m

close all;
clear all;
clc;

% extract time and position
[times,xpos,ypos]=textread('2Ddata1.txt','%s %d %d'); 


[Y, M, D, H, MN, S] = datevec(times);
%convert time string to second
outtime = H*3600+MN*60+S; 
%obtain delta t list
dt_list = diff(outtime);

%obtain delta x and delta y lists
dx = diff(xpos);
dy = diff(ypos);
%obtain velocity lists
xvel = dx./dt_list;
yvel = dy./dt_list;

%store into a measurements matrix
%note that we ignore the first data point (reserve for initialization)
measurements=[xpos(2:end),ypos(2:end), xvel, yvel];

% initializations
w = 0.001;
varx = 1;
vary = 1;

%4x1
x=[xpos(2);ypos(2);xvel(1);yvel(1)];
%The State Transition Matrix phi 4x4
%initialize using the first delta t
phi=[1 0 dt_list(1) 0;
    0 1 0 dt_list(1);
    0 0 1 0;
    0 0 0 1]; 

%Measurement Matrix 2x4
H=[1 0 0 0;0 1 0 0]; 

P =eye(size(phi));
R=[varx 0;0 vary];
%Initial Process Noise Covariance 4x4
%initialize using the first delta t
Q=[(dt_list(1)^4)/4 0 (dt_list(1)^3)/2 0;
    0 (dt_list(1)^4)/4 0 (dt_list(1)^3)/2;
    (dt_list(1)^3)/2 0 (dt_list(1)^2) 0;
    0 (dt_list(1)^3)/2 0 (dt_list(1)^2)]*w; 

Pm = eye(size(phi));
xm = x;

%Initialize container
output = [];

for i=1:length(dt_list)

    %Updating matrix phi with current delta t
    phi=[1 0 dt_list(i) 0;
        0 1 0 dt_list(i);
        0 0 1 0;
        0 0 0 1];
    %Updating Q matrix with current delta t
    Q=[(dt_list(i)^4)/4 0 (dt_list(i)^3)/2 0;
        0 (dt_list(i)^4)/4 0 (dt_list(i)^3)/2;
        (dt_list(i)^3)/2 0 (dt_list(i)^2) 0;
        0 (dt_list(i)^3)/2 0 (dt_list(i)^2)]*w;
        
    %Get current measurement
    z=H*measurements(i,:)';
    %Perform measurement update
    [x, P] = measurement_update(z, H, xm, Pm, R); 
    %Store prediction
    output = [output x];
    
    %Perform time update
    [xm, Pm] = time_update(x, phi, P, Q); 
end

% Display and plots x-y
figure(1);
plot(measurements(:,1),measurements(:,2));
title('2D position (x-y coordinate)');
xlabel('x position (m)');
ylabel('y position (m)');
hold on;
plot(output(1,:),output(2,:));
hold off;

figure(2);
plot(measurements(:,3),measurements(:,4));
title('2D velocity (x-y coordinate)');
xlabel('x velocity (m/s)');
ylabel('y velocity (m/s)');
hold on;
plot(output(3,:),output(4,:));
hold off;
 
% Position vs. time
figure(3);
plot(outtime(2:end), measurements(:,1));
title('x position vs. time');
xlabel('time (s)');
ylabel('x position (m)');
hold on;
plot(outtime(2:end), output(1,:));
hold off;

figure(4);
plot(outtime(2:end), measurements(:,2));
title('y position vs. time');
xlabel('time (s)');
ylabel('y position (m)');
hold on;
plot(outtime(2:end), output(2,:));
hold off;

% Velocity vs. time
figure(5);
scatter(outtime(2:end), measurements(:,3))
title('x velocity vs. time');
xlabel('time (s)');
ylabel('x velocity (m/s)');
hold on;
plot(outtime(2:end), output(3,:));
hold off;

figure(6);
scatter(outtime(2:end), measurements(:,4))
title('y velocity vs. time');
xlabel('time (s)');
ylabel('y velocity (m/s)');
hold on;
plot(outtime(2:end), output(4,:));
hold off;


