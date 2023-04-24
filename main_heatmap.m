run('new_env.m');

%% PROGRAM CONFIGURATION
graph = 0; % To graph the wave equations
save_graph = 0;
save_results = 0; 
save_colormap = 1;
results_filename = '';

%% 1D Meshing
xstart = 0;             % Start pdisoint
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

MAP = zeros(1,50);
map_index = 1;

for E0 = 1 : 1.6 : 9
    for l0 = 0.2 : 0.96 : 5
        for x0 = 3 : 0.8 : 7
            Evble =@(x)E+E0*(heaviside(x-(x0-(l0/2)))-heaviside(x-(x0+(l0/2))));
            run('WaveEquation_ext.m');
            map_index = map_index + 1;
        end
    end
end

colormap('parula')
imagesc(MAP)
colorbar
