%% Some initial stuff
clear, clc, close all
%% Constants
Lx = 10;   % Domain length
Ly = 10;
dx = 0.1; % grid granularity
dy = dx;
nx = fix(Lx/dx); % number of points
ny = fix(Ly/dy);
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);

T = 30; % total time

CFL = 0.5;
c = 1;
dt = CFL*dx/c;

K = 0.1; % dampening constant
%% Domain
h = zeros(nx, ny); % represents field
h_prev = h;
h_next = h;

%% Initial condition
% h(1:3) = [0.1 0.2 0.1];
%% Actual simulation
t=0;
while(t < T)
    % Boundary condition (Reflecting boundary condition)
%     h(:,[1 ny])=0;
%     h([1 nx],:)=0;
    
%     % Absorbing boundary condition
    h_next(1,:)=h(2,:) + (CFL-1)/(CFL+1)*(h_next(2,:)-h(1,:));
    h_next(nx,:)=h(nx-1,:) + (CFL-1)/(CFL+1)*(h_next(nx-1,:)-h(nx,:));
    h_next(:,1)=h(:,2) + (CFL-1)/(CFL+1)*(h_next(:,2)-h(:,1));
    h_next(:,ny)=h(:,ny-1) + (CFL-1)/(CFL+1)*(h_next(:,ny-1)-h(:,ny));
    
    t = t+dt; % stepping time
    h_prev = h;
    h = h_next;
    
    % Source
    if(t < 10)
        h(50,50) = dt^2*20*sin(30*pi*t/20);
    end
    
    for i=2:nx-1
        for j=2:ny-1
            % without dampening
            h_next(i,j) = 1/(1+K*dt)*(2*h(i,j) - h_prev(i,j) + K*dt*h_prev(i,j) + ...
                CFL^2*(h(i+1,j)+h(i,j+1)-4*h(i,j)+h(i-1,j)+h(i,j-1)));
        
%         % with dampening
%         h_next(i) = 1/(1+K*dt)*(2*h(i)-h_prev(i) + K*dt*h_prev(i) + CFL^2*(h(i+1)-2*h(i)+h(i-1)));
        end
    end

    surf(x,y,h,'EdgeColor', 'none');
    axis([0 Lx 0 Ly -0.05 0.05]);
    pause(0.01);
    
end