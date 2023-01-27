%% Visualization for the coupled run
%% Load data
f = fopen( 'solution-coarse.txt', 'r' );
t = 0;
solution_coarse = {};
while ~feof(f)
    line = fgets(f);
    if length(line) >= 8 && strcmp( line(1:8), 'TIMESTEP' )
        t = sscanf( line, '%*s%d' );
        solution_coarse{t+1} = [];
    else
        s = sscanf( line, '%f%f' );
        solution_coarse{t+1} = [ solution_coarse{t+1}; s' ];
    end
end
fclose(f);
f = fopen( 'solution-fine.txt', 'r' );
t = 0;
solution_fine = {};
while ~feof(f)
    line = fgets(f);
    if length(line) >= 8 && strcmp( line(1:8), 'TIMESTEP' )
        t = sscanf( line, '%*s%d' );
        solution_fine{t+1} = [];
    else
        s = sscanf( line, '%f%f' );
        solution_fine{t+1} = [ solution_fine{t+1}; s' ];
    end
end
fclose(f);
for i = 1 : length(solution_fine)
    solution{i} = sortrows( [ solution_fine{i}; solution_coarse{i} ] );
end

%% Create UI
close all
f = figure;
ax = axes( 'Parent', f, 'position', [0.12, 0.2, 0.8, 0.75] );
hold on;

%% Set figure updater callback
b = uicontrol( 'Parent', f, 'Style', 'slider', 'Position', [20,20,500,20], 'value', 0, 'min', 0, 'max', length(solution)-1 );
clear    = @(t) cla;
redraw   = @(t) plot( solution{t}(:,1), solution{t}(:,2) );
scatter1 = @(t) scatter( solution_coarse{t}(:,1), solution_coarse{t}(:,2), 'xr' );
scatter2 = @(t) scatter( solution_fine{t}(:,1), solution_fine{t}(:,2), 'xb' );
relim    = @(t) axis([0,11,0,4]);
b.Callback = @(src,evt) cellfun(@feval, {clear, redraw, scatter1, scatter2, relim}, {0, floor(src.Value) + 1, floor(src.Value) + 1, floor(src.Value) + 1, 0})
cellfun(@feval, {redraw, scatter1, scatter2, relim}, {1, 1, 1, 0})

