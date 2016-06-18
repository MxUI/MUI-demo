%% Visualization for standalone run
%% Load data
f = fopen( 'solution.txt', 'r' );
t = 0;
solution = {};
while ~feof(f)
    line = fgets(f);
    if length(line) >= 8 && strcmp( line(1:8), 'TIMESTEP' )
        t = sscanf( line, '%*s%d' );
        solution{t+1} = [];
    else
        s = sscanf( line, '%f%f' );
        solution{t+1} = [ solution{t+1}; s' ];
    end
end
fclose(f);

%% Create UI
close all
f = figure;
ax = axes( 'Parent', f, 'position', [0.12, 0.2, 0.8, 0.75] );
plot( solution{1}(:,1), solution{1}(:,2) );
ylim( [0,1] );

%% Set figure updater callback
b = uicontrol( 'Parent', f, 'Style', 'slider', 'Position', [20,20,500,20], 'value', 0, 'min', 0, 'max', length(solution) );
redraw = @(t) plot( solution{t}(:,1), solution{t}(:,2) );
relim  = @(t) ylim([0,1]);
b.Callback = @(src,evt) cellfun(@feval, {redraw, relim}, {floor(src.Value) + 1, 0})
