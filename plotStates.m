function [Nstates] = plotStates(states, obs)
    T = length(states);
    assert(T == length(obs));
    
    uniqueStates = unique(states);
    Nstates = length(uniqueStates);
    
    % Number down
    lowStates(Nstates) = 0;
    for i = 1:Nstates
        k = uniqueStates(i);
        lowStates(states == k) = i;
    end

    colors = prism(Nstates);    
    %colors = copper(Nstates);    


    figure;
    hold on;
    
    stateStrs = {};
    for t = 1:T
        scatter(t, obs(t), 8, colors(lowStates(t),:), 'filled');        
    end
        
end
