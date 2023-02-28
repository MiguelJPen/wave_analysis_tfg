run('new_env.m');

%% PROGRAM CONFIGURATION
graph = 0; % To graph the wave equations
save_graph = 1;
save_results = 0; 
results_filename = 'results.csv';

%% 1D Meshing
xstart = 0;             % Start point
xend   = 10;            %%% End point
tne    = 100;           % Total number of element in the domain.
dx     = 10/100;         % spatial step
%%% spatial step dx = 10/100 = 0.1

% Element type:         Q1 --> LINEAR,  Q2 --> QUADRATIC
elementtype = 'Q2';
% Creating 1D Mesh.
[ L, lnn, nne, el, egnn, tnn, x ] = CreateMesh( elementtype, tne, xstart, xend  );
        
%% Material Properties (Constant with in elements -- Q0)

% MECHANICAL
%%% We nondimensionalize
E     = 1;     % Outer velocity ; % Elasticity Tensor
rho   = 1;     % Density of the material

x0    = 5;  % range [2,8]
l0    = 2;  % range [0.5 5.5]
E0    = 9;  % range [1,9]
Evble =@(x)E+E0*(heaviside(x-(x0-l0))-heaviside(x-(x0+l0)));

run('WaveEquation_ext.m');

