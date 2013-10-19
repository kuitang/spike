function plotDeviations(ss, zz, obs)
% plotDeviations(ss, zz, obs) plot states of zz and ss that deviate.
    
    T = length(obs);
    
    [ tab, loTab, cost, dist ] = matchStates(ss, zz);
    
    devRows = tab(:,3) ~= 0;
    devTab  = tab(devRows,:);
    
    devIdxs      = zeros(1, T);
    disIdxs = zeros(1, T);
    
    nDevs = sum(devRows);
    for n = 1:nDevs
        s = devTab(n,1);
        z = devTab(n,2);
        
        sidxs = ss == s;
        zidxs = zz == z;
        
        devIdxs(sidxs) = 1;
        devIdxs(zidxs) = 1;
        
        disIdxs(sidxs ~= zidxs) = 1;
    end
    
    % Wash out the non-deviated states
    devObs = obs;
    devObs(~devIdxs) = NaN;
    
    devSS = ss;
    devSS(~devIdxs) = 0;
    
    devZZ = zz;
    devZZ(~devIdxs) = 0;
    
%     ns = plotStates(devSS, devObs);    
%     plotStates(devZZ, devObs);
    
    disObs = obs;
    disObs(~disIdxs) = NaN;
    
    disSS = ss;
    disSS(~disIdxs) = 0;
    
    disZZ = zz;
    disZZ(~disIdxs) = 0;
    
    nsS = plotStates(disSS, disObs);
    title(['Points with unmatched states; HMM states; ' num2str(nsS) ' states']);
    nsZ = plotStates(disZZ, disObs);    
    title(['Points with unmatched states; DP emission states; ' num2str(nsZ) ' states']);
end

