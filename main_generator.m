run('new_env.m');

%% PROGRAM CONFIGURATION
graph = 0; % To graph the wave equations
save_graph = 0;
save_results = 1; 
save_colormap = 0;
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

if save_results ~= 0
    resFile = fopen(results_filename, 'a+');
end

% MECHANICAL
%%% We nondimensionalize
E     = 1;     % Outer velocity ; %2e5 Elasticity Tensor
rho   = 1;     % 11.6e2 ; % Density of the material

% x0 range [3,   7]
% l0 range [0.2, 5]
% E0 range [1,   9]

for x0 = 3 : 0.4444444 : 7
    for l0 = 0.2 : 0.5333333 : 5
        for E0 = 1 : 0.8888888 : 9
%             if (x0-(l0/2)) >= 2.5 
%                 if (x0+(l0/2)) <= 7.5
                     Evble =@(x)E+E0*(heaviside(x-(x0-(l0/2)))-heaviside(x-(x0+(l0/2))));
                     run('WaveEquation_ext.m');
%                 end
%             end
         end
    end
end

if save_results ~= 0
    fclose(resFile);
end