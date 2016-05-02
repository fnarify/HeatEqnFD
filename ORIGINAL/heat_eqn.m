% nx - number of points in the x-direction
% dt - time step
% nt - number of time steps
%
% Stuart C. Hawkins - 18 May 2015

function heat_eqn(nx,dt,nt)

% set mesh spacings based on nx and ny
dx=1/(nx+1);

% if dt is zero then choose a time step that satisfies the CFL condition
if dt==0    
    dt=0.4*dx^2;  
end

% get mesh points
x=dx*(0:nx+1);
t=dt*(0:nt);

% get interior points... for Dirichlet problem we only work with the
% interior points
xi=x(2:end-1);

% get the parameter in the FD discretisation
nu=dt/dx^2;

% setup the FD matrix
A=diag((1-2*nu)*ones(nx,1))+diag(nu*ones(nx-1,1),1)+diag(nu*ones(nx-1,1),-1);

% initialise the values of the solution... note: this basically reserves
% space for the values and makes the code run faster.
% The columns hold the values at each time step
u_approx=zeros(nx,nt+1);

% put the initial condition into the first column of the matrix
u_approx(:,1)=F0(xi);

% loop through the time steps
for k=1:nt
    % compute the value at the next time step by multiplying by the FD
    % matrix
    u_approx(:,k+1)=A*u_approx(:,k);
end

% visualise the results... the 4th parameter picks out some timesteps to
% plot the solution at
visualise(xi,t,u_approx,[1 round(length(t)/4) round(3*length(t)/4) length(t)],'')

%-----------------------------------------------------
% subfunction that specifies the initial condition
%-----------------------------------------------------

function z=F0(x)

z=(x>=0.25) & (x<=0.75);
