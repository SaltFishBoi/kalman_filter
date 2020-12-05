% main3D.m

close all;
clear all;
clc;

% extract time and position
[times,xpos,ypos,zpos]=textread('3Ddata1.txt','%s %f %f %f');

[Y, M, D, H, MN, S] = datevec(times);
%convert time string to second
outtime = H*3600+MN*60+S;
%obtain delta t list
dt_list = diff(outtime);

%obtain delta x and delta y lists
dx = diff(xpos);
dy = diff(ypos);
dz = diff(zpos);
%obtain velocity lists
xvel = dx./dt_list;
yvel = dy./dt_list;
zvel = dz./dt_list;

%store into a measurements matrix
%note that we ignore the first data point (reserve for initialization)
measurements=[xpos(2:end),ypos(2:end), zpos(2:end), xvel, yvel, zvel];

% initializations
w = 0.0001;
varx = 1000;
vary = 1000;
varz = 1000;

%6x1
x=[xpos(2);ypos(2);zpos(2);xvel(1);yvel(1);zvel(1)];
%The State Transition Matrix phi 6x6
%initialize using the first delta t
phi=[1 0 0 dt_list(1) 0 0;
    0 1 0 0 dt_list(1) 0;
    0 0 1 0 0 dt_list(1);
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1]; 

%Measurement Matrix 3x6
H=[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0]; 

P =eye(size(phi));
R=[varx 0 0; 0 vary 0; 0 0 varz];
%Initial Process Noise Covariance 6x6
%initialize using the first delta t
Q=[(dt_list(1)^4)/4 0 0 (dt_list(1)^3)/2 0 0;
    0 (dt_list(1)^4)/4 0 0 (dt_list(1)^3)/2 0;
    0 0 (dt_list(1)^4)/4 0 0 (dt_list(1)^3)/2;
    (dt_list(1)^3)/2 0 0 (dt_list(1)^2) 0 0;
    0 (dt_list(1)^3)/2 0 0 (dt_list(1)^2) 0;
    0 0 (dt_list(1)^3)/2 0 0 (dt_list(1)^2)]*w;

Pm = eye(size(phi,1));
xm = x;

%Initialize container
output = [];

for i=1:length(dt_list)

    %Updating matrix phi with current delta t
    phi=[1 0 0 dt_list(1) 0 0;
        0 1 0 0 dt_list(1) 0;
        0 0 1 0 0 dt_list(1);
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1];
    %Updating Q matrix with current delta t
    Q=[(dt_list(1)^4)/4 0 0 (dt_list(1)^3)/2 0 0;
        0 (dt_list(1)^4)/4 0 0 (dt_list(1)^3)/2 0;
        0 0 (dt_list(1)^4)/4 0 0 (dt_list(1)^3)/2;
        (dt_list(1)^3)/2 0 0 (dt_list(1)^2) 0 0;
        0 (dt_list(1)^3)/2 0 0 (dt_list(1)^2) 0;
        0 0 (dt_list(1)^3)/2 0 0 (dt_list(1)^2)]*w;
        
    %Get current measurement
    z=H*measurements(i,:)';
    %Perform measurement update
    [x, P] = measurement_update(z, H, xm, Pm, R); 
    %Store prediction
    output = [output x];
    
    [xm, Pm] = time_update(x, phi, P, Q); 
end

% Display and plots
figure(1);
title('3D position (x-y-z coordinate)')
xlabel('x position (yd)');
ylabel('y position (yd)');
zlabel('z position (yd)');
scatter3(measurements(:,1),measurements(:,2), measurements(:,3))
hold on
plot3(output(1,:),output(2,:),output(3,:));
hold off

figure(2);
title('3D position (x-y-z coordinate)')
xlabel('x velocity (yd/s)');
ylabel('y velocity (yd/s)');
zlabel('z velocity (yd/s)');
scatter3(measurements(:,4),measurements(:,5), measurements(:,6))
hold on
plot3(output(4,:),output(5,:),output(6,:));
hold off

% Position vs. time
figure(3);
plot(outtime(2:end), measurements(:,1));
title('x position vs. time');
xlabel('time (s)');
ylabel('x position (yd)');
hold on
plot(outtime(2:end), output(1,:));
hold off

figure(4);
plot(outtime(2:end), measurements(:,2));
title('y position vs. time');
xlabel('time (s)');
ylabel('y position (yd)');
hold on
plot(outtime(2:end), output(2,:));
hold off

figure(5);
plot(outtime(2:end), measurements(:,3));
title('z position vs. time');
xlabel('time (s)');
ylabel('z position (yd)');
hold on
plot(outtime(2:end), output(3,:));
hold off

% Velocity vs. time
figure(6);
scatter(outtime(2:end), measurements(:,4))
title('x velocity vs. time');
xlabel('time (s)');
ylabel('x velocity (yd/s)');
hold on;
plot(outtime(2:end), output(4,:));
hold off;

figure(7);
scatter(outtime(2:end), measurements(:,5))
title('y velocity vs. time');
xlabel('time (s)');
ylabel('y velocity (yd/s)');
hold on;
plot(outtime(2:end), output(5,:));
hold off;

figure(8);
scatter(outtime(2:end), measurements(:,6))
title('z velocity vs. time');
xlabel('time (s)');
ylabel('z velocity (yd/s)');
hold on;
plot(outtime(2:end), output(6,:));
hold off;

