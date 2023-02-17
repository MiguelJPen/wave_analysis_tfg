run('new_env.m');

%% PROGRAM CONFIGURATION
graph = 0; % To graph the wave equations
save_results = 1; 
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
E     = 1;     % Outer velocity ; %2e5 Elasticity Tensor
rho   = 1;     % 11.6e2 ; % Density of the material

%%% Inclusion  30*25*40=30000 parameter sets
% x0 range [2,8] step 0.2 6*5=30 values
% l0 range [0.5 5.5] step 0.2 5*5=25 values
% E0 range [1,9] step 0.2 5*8=40  values

for x0 = 2 : 0.2 : 8
    for l0 = 0.5 : 0.2 : 5.5
        for E0 = 1 : 0.2 : 9
            Evble =@(x)E+E0*(heaviside(x-(x0-l0))-heaviside(x-(x0+l0)));
            run('WaveEquation_ext.m');
        end
    end
end

