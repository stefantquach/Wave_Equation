%% Some initial stuff
clear, clc, close all
%% Constants
L = 10;   % Domain length
dx = 0.05; % grid granularity
nx = fix(L/dx); % number of points
x = linspace(0, L, nx);

T = 30; % total time

CFL = 1;
c = 1;
dt = CFL*dx/c;

K = 0.1; % dampening constant
%% Domain
h = zeros(nx, 1); % represents field
h_prev = h;
h_next = h;

%% Initial condition
% h(1:3) = [0.1 0.2 0.1];
%% Actual simulation
t=0;
while(t < T)
    % Boundary condition (Reflecting boundary condition)
    h([1 nx])=0;
    
%     % Absorbing boundary condition
%     h_next(1)=h(2) + (CFL-1)/(CFL+1)*(h_next(2)-h(1));
%     h_next(nx)=h(nx-1) + (CFL-1)/(CFL+1)*(h_next(nx-1)-h(nx));
    
    t = t+dt; % stepping time
    h_prev = h;
    h = h_next;
    
    % Source
    if(t < 10)
        h(1) = dt^2*30*sin(30*pi*t/T);
    end
    
    for i=2:nx-1
        % without dampening
        % h_next(i) = 2*h(i)- h_prev(i) + CFL^2*(h(i+1)-2*h(i)+h(i-1));
        
        % with dampening
        h_next(i) = 1/(1+K*dt)*(2*h(i)-h_prev(i) + K*dt*h_prev(i) + CFL^2*(h(i+1)-2*h(i)+h(i-1)));
    end

    plot(x, h);
    xlim([0, L]);
    ylim([-1, 1]);
    pause(0.01);
    
end