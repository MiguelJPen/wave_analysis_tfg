%% 1-D ELASTO-DYNAMICS CODE
% 
% Equations: \nabla\cdot\sigma + \rho a = 0
% 	     \sigma: stess, rho: density, a: acceleration
%
% Author: Abdullah Waseem       
% Created: 2019
% Contact: engineerabdullah@ymail.com

%% Pre-calculation of Gauss-Legendre Quadrature, Shape function and their Derivatives
% Gauss Quadrature
nnp = 3;
run('GaussianLegendre.m');
% Shape Functions
run('ShapeFunctions.m');

%% 1D FEM CORE

% Initializing Element Matrices
Me = zeros(nne, nne, tne);      % Mass
Ke = zeros(nne, nne, tne);      % Stiffness
Fe = zeros(nne, 1  , tne);      % Force

% Element loop
for en = 1 : tne
    % Gauss integration loop
	for gs = 1 : ngp
		
		% Jacobian Matrix
		Jcbn = B(gs,:)*x(egnn(en,:));
		% Iso-parameteric map
		x_z  = N(gs,:) * x(egnn(en,:));

		%Force at that gauss point
		%%%force = (3*x_z + x_z^2)*exp(x_z);  % This is an example
        f_amp = 1;
        mu = 1; % the center must be inside the interval
        sigma = 0.2;
        force = f_amp * exp(-((x_z-mu)/sigma)^2);

        %%%density at that gauss point
        density = rho;

        %%%velocity at that gauss point
        %%%velocity = sqrt(E/rho)
        %%% elastic constant at that point
        elastic = Evble(x_z);  % velocity^2
		
		% Element Mass Matrix
		Me(:,:,en) = Me(:,:,en) + N(gs,:)' * density * N(gs,:) * glw(gs) * Jcbn;
        % Element Stiffness Matrix
		Ke(:,:,en) = Ke(:,:,en) + B(gs,:)'/Jcbn * elastic * B(gs,:)/Jcbn * glw(gs) * Jcbn;
		% Element Force Vector
		Fe(:,1,en) = Fe(:,1,en) + N(gs,:)' * force * glw(gs) * Jcbn;
        
	end
end

% Assemble barK, barC and barF
[ barM, barK, barF ] = Assembler( egnn, nne, tne, tnn, Me, Ke, Fe, 'sparse' );

%% BOUNDARY CONDITIONS

% The node-1 is the Dirichlet boundary and the last-node is the Neumann boundary.

%  MECHANICAL/THERMAL -- FREE AT BOTH ENDS
p = [];                  %%% Prescribed
f = setdiff(1:tnn,p);    % Free

% Partitioning the matrices
Kpp = barK(p,p); Kpf = barK(p,f); Kfp = barK(f,p); Kff = barK(f,f);
Mpp = barM(p,p); Mpf = barM(p,f); Mfp = barM(f,p); Mff = barM(f,f);

%% Newmark Family -- A Method

% TIME DATA
T    = 10;                  %%% Total Time
cmax = sqrt(9);             %%% maximum expected velocity
cfl  = 1/2/cmax;            %%% cfl stability condition
dt   = cfl*dx;              %%% Time Step Size
tnts = floor(T/dt+1);       %%% Total Number of Time Steps
titrt= 0 : dt : T;

% Newmark parameters
gamma = 0.5;
beta  = 0.5;

% Frequency dependent function for varying a value in time. 
%%%freq = 1;
%%%a = freq*(2*pi*titrt/T);

% Force vector changing in time.
FV = zeros(tnn, tnts);
Val = 1;           % The magnitude of the applied force.
freq = 1;            % dimensionless frequency.
%%%FV(end,1:floor(tnts/2)) = Val*sin(a(1:floor(tnts/2))); end force
vect = Val*(1-2*pi^2*freq^2*titrt.^2).*exp(-pi^2*freq^2*titrt.^2);
FV(1,1:end) = vect;

% Initializing displacements in time.
U = zeros(tnn, tnts);
% Velocities
V = zeros(tnn, tnts);
% Accelerations.
A = zeros(tnn, tnts);

% Calculating values at time 0 - step 1
A(f,1) = Mff \ (FV(f,1) - Kff*U(f,1));

% When the conductivity and the capacity matrices does not change i.e. linear case
% The system matrix can be assembled, combined and decomposed for faster simulations
% This feature was first introduced in Matlab 2017b. If you have an older version of 
% Matlab then remove the word "decomposition" from the following line.
dK = decomposition(Mff + Kff*dt^2*beta);

if save_results ~= 0
    resFile = fopen(results_filename, 'a+');
    fprintf(resFile, '\n %6d, ', x0);
    fprintf(resFile, '%2.4f, ', l0);
    fprintf(resFile, '%2.4f, ', E0);
end

% Time marching (Newmark Family -- A Method)
for n = 2 : tnts
	
	% Displacement Predictor
	U(f,n) = U(f,n-1) + dt*V(f,n-1) + dt^2/2*(1-2*beta)*A(f,n-1);
    
	% Velocity Predictor
	V(f,n) = V(f,n-1) + dt*(1-gamma)*A(f,n-1);
	
	% Solution
	A(f,n) = dK \ (FV(f,n) - Kff*U(f,n));
	
	% Displacement Corrector
	U(f,n) = U(f,n) + dt^2*beta*A(f,n);
	% Velocity Corrector
	V(f,n) = V(f,n) + dt*gamma*A(f,n);

    if save_results ~= 0
        fprintf(resFile, '%2.5f, ', A(tnn, n));
    end
	
    if graph ~= 0
        figure(1);
        subplot(3,1,1)
        plot(x,U(:,n));
        xlim([min(x)-max(x)/10 max(x)+max(x)/10])
        ylim([-1 1]/5e0)
        title('Displacement');
        drawnow

        subplot(3,1,2);
        plot(x,V(:,n));
        xlim([min(x)-max(x)/10 max(x)+max(x)/10])
        ylim([-1 1]/1e0)
        title('Velocity');
        drawnow

        subplot(3,1,3)
        plot(x,A(:,n));
        xlim([min(x)-max(x)/10 max(x)+max(x)/10])
        ylim([-1 1]/1e-1)
        title('Acceleration');
        drawnow
    end
end

if save_results ~= 0
    fclose(resFile);
end