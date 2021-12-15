function [xref,uref] = computeReference(theta,tend,N,period)
%computeReference computes the periodic persistently exciting reference
%trajectory around a desired steady state given the uncertain parameter
%theta.
% theta: uncertain parameter
% tend: MPC finite horizon
% N: simulation time
% period: time period for reference 
% 

% Define sufficient desired steady-state values
% xs = [xs1 xs2 ... xsn]
% us = [us1 us2 ... usm]
% and define the rest as symbolic variables such that the related system 
% of equations 
%   xs = f(xs,us,theta)
% can be solved via the 'solve' command.


% Generate trajectory
ubar = ;  % ampltiude around steady us

% Generate sinosoids as persistently exciting control sequence
% Adjust to number of inputs
for k=1:period
    uper(k)=us+ubar*sin(2*pi/period*(k-1));
end

%% Compute related periodic state sequence
syms x1 x2
x=[x1;x2];

% Recursive dynamics
for k=1:period
    x=dynamics(x, uper(k), 0, theta);
end

% Solve for x_r(0)
[x1,x2] = vpasolve(x(1)-x1,x(2)-x2,x1,x2,xs);

% Generate x sequence
xper=[double(x1);double(x2)];
for k=1:period-1
    xper(:,k+1)=x_dyn_reactor(xper(:,k), uper(k), 0, theta);
end

%% Extend x,u until tend
nperiods   = floor(tend/period)+1+N;
xref       = repmat(xper,1,nperiods);
uref       = repmat(uper,1,nperiods);