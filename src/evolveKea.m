function out = evolveKea( keavars,inputs,params )

OSM_init = params.LumenClStart + params.LumenKStart + params.LumenMgStart;
Volume_init = params.LumenVolumeStart;

LumenMg = inputs.LumenMg;
LumenCl = inputs.LumenCl;
LumenProtons = keavars(1,:);
LumenK = keavars(2,:);
Pkea = params.Pkea;

volume = keavars(3,:);

lumen.Protons = keavars(1,:);
lumen.Mg = inputs.LumenMg;
lumen.Cl = inputs.LumenCl;
lumen.K = keavars(2,:);
stroma.Protons=params.StromaProtonsStart;
stroma.Cl=params.StromaClStart;
stroma.Mg=params.StromaMgStart;
stroma.K=params.StromaKStart;

s=getStaticThylakoidValues(lumen, stroma, params);



rates(1,:)=getflux(lumen.Mg,  params.StromaMgStart,  s.deltamuMg,s.deltapsi,  params.PMg, params.zMg);
rates(2,:)=getflux(LumenCl,  params.StromaClStart,  s.deltamuCl,s.deltapsi,  params.PCl, params.zCl);
rates(3,:)=getflux(LumenK,   params.StromaKStart,   s.deltamuK, s.deltapsi,  params.PK, params.zK);

dLumenMg        =  fluxToConcentration(rates(1,:), 'lumen', params,volume);
dLumenCl        =  fluxToConcentration(rates(2,:), 'lumen', params,volume);
dLumenK         =  fluxToConcentration(rates(3,:), 'lumen', params,volume);


rate = ( -LumenProtons .* stroma.K .* Pkea ) + ( stroma.Protons .* dLumenK .* Pkea );
dLumenProtons = rate;
dLumenK = -rate;

OSM = getOsmolarity(dLumenK,dLumenCl,dLumenMg);

Volume = ( OSM .* Volume_init ) ./ OSM_init;

out = [ dLumenProtons;dLumenK;Volume];

%% Load in the variables we'll be changing
% LumenProtons = keavars( 1,: );
% LumenK = keavars( 2,: );
% 
% 
% 
% %% Some static variables
% StromaK = params.StromaKStart;
% StromaProtons = params.StromaProtonsStart;
% Pkea = params.Pkea;
% 
% rate = (-LumenProtons .* StromaK .* Pkea) + (StromaProtons .* LumenK .* Pkea);
% 
% dLumenProtons = rate;
% dLumenK = -rate;
% 
% 
% out = [ dLumenProtons;dLumenK];













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

