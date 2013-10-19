function [ tab, loTab, cost, dist ] = matchStates(x, y)
% [tab, leftX, leftY] = matchStates(x, y) find closest isomorphism of state sequences.
%
%   Assume each state has unique count for now; we eventually need to relax
%   to something more discriminatory, but note that generally isomorphism
%   would be NP-hard.
%
%   Requires rectangular assignment problem toolbox; uses Munkres'
%   algorithm.

    [dx, xLoToHi, Nx] = numberDown(x);
    [dy, yLoToHi, Ny] = numberDown(y);
    
    % Cost is the Hamming distance of the indicator vectors for each state.
    % We use the down-numbered states, so states are consecutive integers.
    dist(Nx,Ny) = 0;
    for nx = 1:Nx        
        for ny = 1:Ny
            nxIdxs = dx == nx;
            nyIdxs = dy == ny;
            
            dist(nx, ny) = sum(nxIdxs ~= nyIdxs);
        end
    end    
    
    %assert(all(vec(dist == dist')), 'dist was not symmetric');
    
    %dist(end+1,end) = 0;

    [assignment, cost] = assignmentoptimal(dist);
    % assignment(x) = y.
    
    tab = assignment;
    
    % Now find the true numbers
    for i = 1:Nx        
        loTab(i,1) = i;
        loTab(i,2) = assignment(i);
        loTab(i,3) = dist(i, assignment(i));
        
        tab(i,1) = xLoToHi(i);
        tab(i,2) = yLoToHi(assignment(i));               
        tab(i,3) = dist(i, assignment(i));
    end
end
