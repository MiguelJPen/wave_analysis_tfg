run('new_env.m');

%% PROGRAM CONFIGURATION
graph = 0; % To graph the wave equations
save_graph = 0;
save_results = 1; 
results_filename = 'results_lhs.csv';

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
E     = 1;     % Outer velocity ; %2e5 Elasticity Tensor
rho   = 1;     % 11.6e2 ; % Density of the material

%%% Latin hypercube sampling
% x0 range [2,8]
% l0 range [0.5 5.5]
% E0 range [1,9]

n_samples = 10;
X = lhsdesign(n_samples, 3);

x0_l = zeros(n_samples, 1);
l0_l = zeros(n_samples, 1);
E0_l = zeros(n_samples, 1);

for i = 1 : n_samples
    x0_l(i) = 2 + X(i, 1) * 6;
    l0_l(i) = 0.5 + X(i, 2) * 5;
    E0_l(i) = 1 + X(i, 3) * 8;
end

for i = 1 : n_samples
    x0 = x0_l(i);
    l0 = l0_l(i);
    E0 = E0_l(i);
    Evble =@(x)E+E0*(heaviside(x-(x0-l0))-heaviside(x-(x0+l0)));
    run('WaveEquation_ext.m');
end

