
function [ outputSim ] = postProcesses( inputSim )
outputSim = inputSim;
vars = inputSim.simparams.varsforsim;
values = inputSim.simulatedvalues;

% Get Lumen values
LumenClPos = ismember(inputSim.simparams.varsforsim,cellstr('LumenCl'));
LumenMgPos = ismember(inputSim.simparams.varsforsim,cellstr('LumenMg'));
LumenKPos = ismember(inputSim.simparams.varsforsim,cellstr('LumenK'));
LumenProtonsPos = ismember(inputSim.simparams.varsforsim,cellstr('LumenProtons'));

LumenCl = values(LumenClPos,:);
LumenMg = values(LumenMgPos,:);
LumenK = values(LumenKPos,:);
LumenProtons = values(LumenProtonsPos,:);

% Get the stroma values
StromaCl = inputSim.params.StromaClStart;
StromaMg = inputSim.params.StromaMgStart;
StromaK = inputSim.params.StromaKStart;

dOSM = (LumenCl+LumenMg+LumenK)-(StromaCl+StromaMg+StromaK);
% dOSM = (LumenCl+LumenMg+LumenK);

vars{end+1} = 'OSM';
values(end+1,:) = dOSM;

outputSim.simparams.varsforsim = vars;
outputSim.simulatedvalues = values;

end

