function dConc=fluxToConcentration(flux, space, params,volume)
if( nargin < 4 )
    volume = params.lumenVolumePerArea;
end

%flux is in units of moles per 
switch(space)
    case 'lumen'
        volumecorrection=volume;
        direction=1;
    case 'stroma'
        volumecorrection=(volume*params.StromaVolume/params.LumenVolume);
        direction=-1;
end
dConc=direction*flux./volumecorrection;
end