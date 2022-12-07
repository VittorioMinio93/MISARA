%%%%%%%%%%%%% Search all seismic traces%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% February 2021 
% Vittorio Minio
%
%%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path: Input folder of the seismic traces
% staz: Names of the stations
% compo: component/channel used (Z,N,E,F)
% date_search: time reference string 
% comp_chan: station system used (component/channel)
% ff2: list of the all seismic traces

function  ff2=searchAllStaz(path,staz,compo,date_search,comp_chan)
    nsta=length(staz);%%number of stations
    ff2=[];
    %%Execute one of two groups of statements (component/channel)
    switch comp_chan
        case 'Comp'
            %%loop through nsta stations
            for kk=1:nsta
                sta=staz{kk};%%kk-th station
                list=dir(strcat(path,sta,'*BH',compo,'*',date_search,'*.mat'));%%search seismic trace
                 ff2=[ff2;list];%%append seismic trace to the ff2 list
            end
        case 'Chan'
            %%loop through nsta stations
            for kk=1:nsta
                sta=staz{kk};%%kk-th station
                list=dir(strcat(path,sta,'*',date_search,'*.mat'));%%search seismic trace
                 ff2=[ff2;list];%%append seismic trace to the ff2 list
            end
    end
end