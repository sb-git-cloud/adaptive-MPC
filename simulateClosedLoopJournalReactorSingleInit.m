% Template for closed-loop simulation
% Author: Sven Brüggemann, UCSD
% December 2021

% Related to: S. Brüggemann, R. R. Bitmead, 'Forward-looking persistent 
% excitation in model predictive control', Automatica, Volume 136, 2022,
% 110033, ISSN 0005-1098, https://doi.org/10.1016/j.automatica.2021.110033.

%% Init
clear all
close all

% Simulation parameters
sim_time    = 100;
options = optimset('display','off');

% System parameters
nx = ; % # states
m = ;  % # inputs
p = ;  % # outputs
npar = ;% # of uncertain parameters

% RLS
lambda = 0.9; % forgetting factor
theta0 = []; % actual value of initial uncertaint (npar x 1) (may vary w/ time)
theta_init = [; ]; % initial estimation (npar x 1)
period = ;  % periodicity of persistently excting trajectory
T = [;];  % Weighting matrix for recursive least squares estimator (npar x npar)

% Optimization parameters for MPC
np          = ;    % Prediction horizon 
Q	        = ;    % Weight states (nx x nx)
R           = ;    % Weighing control (m x m)

% Process noise (nx)
wmax = [;];
w    = [-wmax(1) + 2*wmax(1)*rand(1,sim_time);
    -wmax(2) + 2*wmax(2)*rand(1,sim_time);
    ...];
wbar = wmax;

% Initial Reference
[x_ref, u_ref] = computeReferenceReactor(theta_init, np, sim_time, period);
x_init      = [];  % initial state
P_1         = [];  % initial covariance (npar x npar)

% Solution initialization
Pk(:,:,1)   = P_1;
u_seq       = zeros(m,np); % Arbitrary initial sequence
t = 0:1:sim_time;

%% Start simulation
for k=1:size(x_init,2)
    xk{k}(:,1)     = x_init(:,k); % Solution vector state x[k]
    xrefk{k}(:,1)  = x_ref(:,1); % Time-varying reference trajectory
    urefk{k}(:,1)  = u_ref(:,1); % Time-varying reference trajectory (u)
    thetak{k}(:,1) = theta_init(:,k); % Solution theta estimate
    tilde_theta{k}(1) = norm(thetak{k}(:,1)-theta0);  % Error norm
    
    fprintf(1,'Time instance i =     ');
    for i=1:sim_time

        fprintf(repmat('\b',1,5));
        fprintf('%4d\n',i);

        % MPC
        [u,fvalms(i),eflag,output] = fmincon(@(u) costfunction(np, xk{k}(:,i), u, ...
             Q, R, x_ref(:,i:i+np-1), u_ref(:,i:i+np-1), thetak{k}(:,end)), ...
             u_seq, [],[],[],[],[],[],[],options);

        % Save first control input
        uk{k}(:,i) = u(1);

        % Use tail for next iteration
        u_seq   = [u(:,2:end) zeros(m,1)];

        % Time-varying parameter
        theta0(:,i+1) = ;
        

        % Apply control to system
        xk{k}(:,i+1)  = dynamics(xk{k}(:,i), uk{k}(:,i), w(:,i), theta0(:,i+1));

        % Recusrive least squares estimator
        % [xbar^+ = x^+ - f0(x,u)]

        xbar(:,i+1)     = xk{k}(:,i+1)-();
        varphi(:,:,i)   = ;

        Dk = lambda*T+varphi(:,:,i)'*Pk(:,:,i)*varphi(:,:,i);
        thetak{k}(:,i+1) = thetak{k}(:,i) ...
                        + Pk(:,:,i)*varphi(:,:,i)/Dk ...
                            *(xbar(:,i+1)-varphi(:,:,i)'*thetak{k}(:,i));

        Pk(:,:,i+1) = inv(lambda*eye(npar)/Pk(:,:,i)+varphi(:,:,i)/T*varphi(:,:,i)');

        tilde_theta{k}(i+1)= norm(thetak{k}(:,i+1)-theta0(:,i+1));

        % Update trajectory
        [x_ref, u_ref] = computeReference(thetak{k}(:,i+1), np, sim_time, period);
        
        xrefk{k}(:,i+1)=x_ref(:,i+1);
        urefk{k}(:,i+1)=u_ref(:,i+1);
    end
    
end


%% Functions
%%%%%%%%%%%%%%%

function cost = costfunction(np, x0, u, Q, R, ...
                x_ref, u_ref, theta)
    cost = 0;
    x = zeros(2, np);
    x = computeOpenloopSolutionReactor(np, x0, u, theta);
    
    for k=1:np
        cost = cost+runningcosts(Q, R, x(:,k)-x_ref(:,k), ...
            u(:,k)-u_ref(:,k));
    end
end

function cost = runningcosts(Q, R, x, u)
    cost = x'*Q*x + u'*R*u;
end

function x = computeOpenloopSolutionReactor(np, x0, u, theta)

    x(:,1) = x0;
    for k=1:np-1
        % Dynamics x
        x(:,k+1) = dynamics(x(:,k), u(k), 0, theta);
    end
end