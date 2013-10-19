function [ lowStates, loToHi, Nstates ] = numberDown(states)

    T = length(states);    
    
    uniqueStates = unique(states);
    Nstates = length(uniqueStates);
    
    % Number down
    lowStates(Nstates) = 0;
    loToHi(Nstates)    = 0;
    for i = 1:Nstates
        k = uniqueStates(i);
        lowStates(states == k) = i;
        loToHi(i) = k;
    end
    
end

