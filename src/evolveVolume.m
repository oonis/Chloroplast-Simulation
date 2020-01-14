function out = evolveVolume( OSMvars,inputs,params )

%% Get the values that we ACTUALLY want
% LumenK = inputs.LumenK;
% LumenCl = inputs.LumenCl;
% LumenMg = inputs.LumenMg;
% LumenProtons = inputs.LumenProtons;
% LumenFreeProtons = getLumenpH( LumenProtons, params );
% LumenFreeProtons = 10.^( -1 * LumenFreeProtons );

OSM_init = params.LumenClStart + params.LumenKStart + params.LumenMgStart;
Volume_init = params.LumenVolumeStart;

LumenK = OSMvars(1,:);
LumenCl = OSMvars(2,:);
LumenMg = OSMvars(3,:);
LumenProtons = OSMvars(4,:);

% LumenFreeProtons = getLumenpH( LumenProtons, params );
% LumenFreeProtons = 10.^( -1 * LumenFreeProtons );

% Stroma values are assumed to be constant
% StromaK = params.StromaKStart;
% StromaCl = params.StromaClStart;
% StromaMg = params.StromaMgStart;

lumen.Protons = LumenProtons;
lumen.Mg = LumenMg;
lumen.Cl = LumenCl;
lumen.K = LumenK;
stroma.Protons=params.StromaProtonsStart;
stroma.Cl=params.StromaClStart;
stroma.Mg=params.StromaMgStart;
stroma.K=params.StromaKStart;

s=getStaticThylakoidValues(lumen, stroma, params);

rates(1,:)=getflux(LumenMg,  params.StromaMgStart,  s.deltamuMg,s.deltapsi,  params.PMg, params.zMg);
rates(2,:)=getflux(LumenCl,  params.StromaClStart,  s.deltamuCl,s.deltapsi,  params.PCl, params.zCl);
rates(3,:)=getflux(LumenK,   params.StromaKStart,   s.deltamuK, s.deltapsi,  params.PK, params.zK);

dLumenMg        =  fluxToConcentration(rates(1,:), 'lumen', params);
dLumenCl        =  fluxToConcentration(rates(2,:), 'lumen', params);
dLumenK         =  fluxToConcentration(rates(3,:), 'lumen', params);

% dStromaMg = fluxToConcentration( StromaMg,'stroma',params );
% dStromaCl = fluxToConcentration( StromaCl,'stroma',params );
% dStromaK = fluxToConcentration( StromaK,'stroma',params );


OSM = getOsmolarity( dLumenK,dLumenCl,dLumenMg );

% dLumenK = [];
dLumenMg = []; dLumenCl = [];
dLumenProtons = [];

dLumenK( 1:numel(LumenK) ) = 0;
dLumenCl( 1:numel(LumenK) ) = 0;
dLumenMg( 1 : numel(LumenK) ) = 0;
dLumenProtons(1:numel(LumenK)) = 0;


% Volume: V_new / V_init = OSM_new / OSM_init  -> V_new = (OSM_new * V_init) / OSM_init
Volume = ( OSM .* Volume_init ) ./ OSM_init;

% out = LumenK;
out = [dLumenK;dLumenCl;dLumenMg;dLumenProtons;Volume];

%% Model For ion Flux
% 
%%

    function FluxLumen=getflux(LumenConc, StromaConc, deltamu,deltapsi, permeability, z)
        
        %provided units: Molar, Molar, Volts, cm/z, none
        %units:    flux out: moles ions/cm^2=moles ion/liter * cm/second *
        %Volts/volts
        %I liter=1000cm^3 so if LumenConc is in
        %moles/liter*1liter/1000cm^3=1e-3moles/cm^3
        flowout=[deltamu]>[0 ];
        %flowout=flows(1)*flows(2); %deltamu and z have the same sign;
        literspercc=1e-3;
        FluxOut=  LumenConc .* literspercc .* permeability    .*          deltamu /params.voltsperlog;
        FluxIn  =  StromaConc .*  literspercc .* permeability  .*        deltamu/params.voltsperlog;
        %
        fluxtype='lin';
        switch fluxtype
            case 'gh'
%                Goldman-Hodgkin-Katz equation
                f=exp(z*deltapsi/params.voltsperlog);
                if (f-1)<1e-10*ones(size(f))
                    FluxLumen=0;
                else
                    
                    FluxLumen=z/params.voltsperlog * permeability *deltapsi.* (StromaConc-LumenConc.*f)./(1-f);
                end
            case 'lin'
                FluxLumen=-(LumenConc.*flowout +StromaConc.*(1-flowout))  .*literspercc.*permeability .* deltamu./params.voltsperlog;
                
        end
    end

end

