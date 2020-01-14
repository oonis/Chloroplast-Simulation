function getAllValues( sim,fileLocation )
% numColumns = numel( sim.lightIntensity ); % This may not even matter

%% Get the static variables

% Plot time first
bigMatrix(1,1) = cellstr('Time');
bigMatrix(1,2:1+numel(sim.timevalues)) = num2cell(sim.timevalues);

[s,~] = getStaticVals( sim );
staticLabels = fields( s );
for x = 1 : numel( staticLabels )
    bigMatrix(end+1,1) = staticLabels(x);
    currValues = s.(staticLabels{x});
    currValues = num2cell( currValues );
    
    bigMatrix(end,2:1+numel(currValues)) = currValues;
end
    
%% Get the dynamic variables

% Plot light intensity seperately
bigMatrix(end+1,1) = cellstr('LightIntensity');
bigMatrix(end,2:1+numel(sim.LightIntensity)) = num2cell( sim.LightIntensity );

dynamicLabels = sim.simparams.varsforsim;
dynamicValues = num2cell( sim.simulatedvalues );

for x = 1 : numel( dynamicLabels )
    bigMatrix(end+1,1) = dynamicLabels(x);
    currValues = dynamicValues(x,:);
    bigMatrix(end,2:1+numel(currValues)) = currValues;
end

%% Print the matrix to a text file

% So I can keep scientific notation
bigMatrix(:,2:end) = cellfun( @num2str,bigMatrix(:,2:end),'UniformOutput',false );

bigMatrix = bigMatrix'; % Transpose for use in Origin

formatString = '%s';
for x = 1 : numel( bigMatrix(1,2:end) )
    formatString = horzcat(formatString,'\t%s');
end
formatString = horzcat(formatString,'\n');


fid = fopen( fileLocation,'wt' );

for x = 1 : numel( bigMatrix(:,1) )
    fprintf( fid,formatString,bigMatrix{x,:} );
end

fclose(fid);
end

