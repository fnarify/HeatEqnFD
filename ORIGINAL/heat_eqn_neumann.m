% nx - number of points in the x-direction
% dt - time step
% nt - number of time steps
%
% Stuart C. Hawkins - 18 May 2015

% Implicit implementation of solution to heat equation.
function heat_eqn_neumann(nx, dt, nt)

% Set mesh spacings based on nx and ny.
dx = 1/(nx + 1);

% If dt is zero then choose a time step that satisfies the CFL condition.
if dt == 0    
    dt = 0.4*dx^2;
end

% Get mesh points.
x = dx*(0: nx + 1);
t = dt*(0: nt);

% Use all points due to Neumann conditions.
xi = x(1: end);

% Get the parameter in the FD discretisation.
nu = dt/dx^2;

% Setup the FD matrix, for an implicit scheme.
A = diag((1 + 2*nu)*ones(nx, 1)) + diag(-nu*ones(nx - 1, 1), 1) + ...
    diag(-nu*ones(nx - 1, 1), -1);
% Add change to matrix given by Neumann boundary conditions. 
A(1, 2) = -2*nu; A(nx, nx - 1) = -2*nu;

% initialise the values of the solution.
% The columns hold the values at each time step.
u_approx = zeros(nx, nt + 1);

% put the initial condition into the first column of the matrix.
u_approx(:, 1) = F0(xi);

% Loop through the time steps.
for k = 1: nt
    % Given system u_m = A*u_m+1, we solve for u_m+1 = A^-1*u_m.
    u_approx(:, k + 1) = A\u_approx(:, k);
end

% Visualise the results. 
% The 4th parameter picks out some timesteps to plot the solution at.
visualise(xi, t, u_approx, [1, round(length(t)/4), round(3*length(t)/4), length(t)], '');

%-----------------------------------------------------
% Subfunction that specifies the initial condition.
%-----------------------------------------------------

function z = F0(x)

z = (x >= 0.25) & (x <= 0.75);
