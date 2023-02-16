%% 1-D ELASTO-DYNAMICS CODE
% 
% Equations: \nabla\cdot\sigma + \rho a = 0
% 	     \sigma: stress
%        \rho: density
%        a: acceleration
%
% Author: Abdullah Waseem       
% Created: 2019
% Contact: engineerabdullah@ymail.com

clear; clc; clf; path(pathdef); format long
addpath FECore/

%% Parameters

% 1D Mesh control
xstart = 0;              % Start point
xend   = 10;             % End point
tne    = 200;            % Total number of element in the domain.

% Element type:         Q1 --> LINEAR,  Q2 --> QUADRATIC
elementtype = 'Q2';

% Material Properties (Constant with in elements -- Q0)
E     = 2e5;      % Elasticity Tensor
rho   = 11.6e2  ;          % Density of the material

% Gauss Quadrature
nnp = 4;

% Force at that gauss point
f_amp = 1e4;
mu = 1;
sigma = 0.2;

% Time simulation data
T    = 2;                 % Total Time
dt   = 0.01;              % Time Step Size
tnts = T/dt+1;            % Total Number of Time Steps
titrt= 0 : dt : T;

% Newmark parameters
% Determines the accuracy of the solution and how it approximates the true 
% solution. If we increase gamma, we increase accuracy, but also
% computation time.
gamma = 0.75;
% Adjust the level of damping in the solution, preventing numerical 
% instability
beta  = 0.75;

%% 1D Meshing

% Creating 1D Mesh.
[ L, lnn, nne, el, egnn, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  );
% x = Nodal Coordinates Vector        
% egnn = Element Global Node Numbering

%% Pre-calculation of Gauss-Legendre Quadrature, Shape function and their Derivatives
run('GaussianLegendre.m');
% Shape Functions
run('ShapeFunctions.m');

%% 1D FEM (Finite Element Method) CORE

% nne = Number of Nodes in the Element
% tne = Total Number of Elements
% Initializing Element Matrices
Me = zeros(nne, nne, tne);      % Mass
Ke = zeros(nne, nne, tne);      % Stiffness
Fe = zeros(nne, 1  , tne);      % Force

% Element loop
for en = 1 : tne % Recorremos los elementos de la malla.
    % Calculamos las matrices Me, Ke y Fe usando la integración de Gauss.
	for gs = 1 : ngp
		
		% Jacobian Matrix
		Jcbn = B(gs,:) * x(egnn(en,:));
		% Iso-parameteric map
		x_z  = N(gs,:) * x(egnn(en,:));

        % Force vector
        force = f_amp * exp(-((x_z-mu)/sigma)^2);
	    
		% Element Mass Matrix
		Me(:,:,en) = Me(:,:,en) + N(gs,:)'      * rho     * N(gs,:)      * glw(gs) * Jcbn;
        % Element Stiffness Matrix
		Ke(:,:,en) = Ke(:,:,en) + B(gs,:)'/Jcbn * E       * B(gs,:)/Jcbn * glw(gs) * Jcbn;
		% Element Force Vector
		Fe(:,1,en) = Fe(:,1,en) + N(gs,:)' * force * glw(gs) * Jcbn;
        
	end
end

% Assemble barK, barC and barF
% barM: global mass matrix
% barK: global stiffness matrix
% barF: global force vector
% La función suma todas las matrices de elementos individuales para formar
% las matrices globales y las devuelve.
[ barM, barK, barF ] = Assembler( egnn, nne, tne, tnn, Me, Ke, Fe, 'sparse' );

%% BOUNDARY CONDITIONS

% The node-1 is the Dirichlet boundary and the last-node is the Neumann boundary.

% Estas condiciones ayudan a ver lo que ocurre en los extremos de la onda,
% ya que la ecuación de onda es una derivada parcial que describe cómo
% viaja la onda por el medio y esto nos ayuda a ver qué ocurre al final.
% Son útiles para modelar el comportamiento y realizar predicciones a lo
% largo del tiempo.

%  MECHANICAL/THERMAL -- FIXED AT BOTH ENDS         
p = 1;                   % Grado de libertad fijo 
f = setdiff(1:tnn,p);    % Conjunto del resto de grados de libertad que no están fijos.

% Partitioning the matrices
% Estas matrices ayudan a determinar las fuerzas y deformaciones en la
% onda. Sirven para describir las interacciones entre los distintos
% elementos de la estructura, como las relaciones entre las fuerzas y
% deformaciones de elementos junto a sus elementos vecinos.
Kpp = barK(p,p); Kpf = barK(p,f); Kfp = barK(f,p); Kff = barK(f,f);
% Representa la matriz de masa en el método Newmark. Esta matriz se usa
% para describir el comportamiento del sistema y cómo responde a cambios
% en la velocidad y la aceleración.
Mpp = barM(p,p); Mpf = barM(p,f); Mfp = barM(f,p); Mff = barM(f,f);

%% Newmark Family -- A Method

% Initializing displacements in time.
U = zeros(tnn, tnts);
% Velocities
V = zeros(tnn, tnts);
% Accelerations.
A = zeros(tnn, tnts);

% Calculating values at time 0 - step 1
A(f,1) = Mff \ (barF(1) - Kff*U(f,1));

% When the conductivity and the capacity matrices does not change i.e. linear case
% The system matrix can be assembled, combined and decomposed for faster simulations
% This feature was first introduced in Matlab 2017b. If you have an older version of 
% Matlab then remove the word "decomposition" from the following line.
dK = decomposition(Mff + Kff*dt^2*beta);

resFile = fopen('results.csv', 'w');
fprintf(resFile, '%6d,\n', f_amp);
fprintf(resFile, '%2.4f,\n', mu);
fprintf(resFile, '%2.4f,\n\n', sigma);

% Time marching (Newmark Family -- A Method)
for n = 2 : tnts
	
	% Displacement Predictor
	U(f,n) = U(f,n-1) + dt*V(f,n-1) + dt^2/2*(1-2*beta)*A(f,n-1);
    
	% Velocity Predictor
	V(f,n) = V(f,n-1) + dt*(1-gamma)*A(f,n-1);
	
	% Solution
	A(f,n) = dK \ (barF(n) - Kff*U(f,n));
	
	% Displacement Corrector
	U(f,n) = U(f,n) + dt^2*beta*A(f,n);
	% Velocity Corrector
	V(f,n) = V(f,n) + dt*gamma*A(f,n);
	
    % POST-PROCESSING
    figure(1);
    subplot(3,1,1)
    plot(x,U(:,n));
    xlim([min(x)-max(x)/10 max(x)+max(x)/10])
    ylim([-1 1])
    title('Displacement');
    drawnow
    subplot(3,1,2);
    plot(x,V(:,n));
    xlim([min(x)-max(x)/10 max(x)+max(x)/10])
    ylim([-1.5 1.5])
    title('Velocity');
    drawnow
    subplot(3,1,3)
    plot(x,A(:,n));
    xlim([min(x)-max(x)/10 max(x)+max(x)/10])
    ylim([-10 10])
    title('Acceleration');
    drawnow
    
    fprintf(resFile, '%3.5f,\n', A(tnn, n));
    
end

fclose(resFile);