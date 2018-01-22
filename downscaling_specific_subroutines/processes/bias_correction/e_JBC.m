function [vecSimOut, vecSimOutV] = e_JBC(vecRefHist, vecRefHistV, vecSimHist, vecSimHistV, vecSim2Bc, vecSim2BcV, type, typeV)
% 
% %
% %
% %
% %FOR TESTING ONLY (REMOVE!!!)
% vecSim2Bc  = vecSimHist;
% vecSim2BcV = vecSimHistV;
% %
% %
% %

matSim2Bc  = [ vecSim2Bc(:),  vecSim2BcV(:)];
matSimHist = [vecSimHist(:), vecSimHistV(:)];
matRefHist = [vecRefHist(:), vecRefHistV(:)];

%Initialize output
vecSimOut = nan(size(vecSim2Bc));
vecSimOutV = vecSimOut;

limit = 3; %Limit for bias correction factor
  
if ~any2d(isnan(matSimHist)) && ~any2d(isnan(matRefHist))
    %Find copulas for training bias correction:
%     ecSimHist = e_cop(matSimHist);
%     ecRefHist = e_cop(matRefHist);
   

%     ecSim2Bc = e_invcop(matSimHist, matSim2Bc);
    mapSim2Ref = e_cop_map(matSimHist, matRefHist);
    mapSimBc2Sim = e_invcop_map(matSim2Bc, matSimHist);

    for kk = 1 : numel(matSim2Bc(:,1))
%         [~, indRefHist] = min(abs(ecRefHist - ecSim2Bc(kk)));
%         [~, indSimHist] = min(abs(ecSimHist - ecSim2Bc(kk)));
        
        indSim2Bc =  mapSimBc2Sim(kk, 1);
        indSimHist = mapSimBc2Sim(kk, 2);
        indRefHist = mapSim2Ref(mapSim2Ref(:, 1) == indSimHist, 2);
        
        if regexpbl(type, 'mult')
            bcFactor = vecRefHist(indRefHist) / vecSimHist(indSimHist);
            bcFactor = max([min([bcFactor, limit]), -limit]);
            
            vecSimOut(indSim2Bc)  =  vecSim2Bc(indSim2Bc) * bcFactor;
        elseif regexpbl(type, 'add')
            vecSimOut(indSim2Bc)  =  vecSim2Bc(indSim2Bc) +  vecRefHist(indRefHist) -  vecSimHist(indSimHist);
        else
            error('eJbc:unknownType',['Type ' type ' has not been programmed for.']);
        end

        if regexpbl(typeV, 'mult')
            bcFactorV = vecRefHistV(indRefHist) / vecSimHistV(indSimHist);
            bcFactorV = max([min([bcFactorV, limit]), -limit]);
            
            vecSimOutV(indSim2Bc) = vecSim2BcV(indSim2Bc) * bcFactorV;
        elseif regexpbl(typeV, 'add')
            vecSimOutV(indSim2Bc) = vecSim2BcV(indSim2Bc) + vecRefHistV(indRefHist) - vecSimHistV(indSimHist);
        else
            error('eJbc:unknownTypeV',['Type ' typeV ' has not been programmed for.']);
        end
        
    end
    clear kk
    
    %For testing:
%     corrIn = corr(matSim2Bc);
%     corrRef = corr(matRefHist);
%     corrOut = corr([vecSimOut,vecSimOutV]);
%     disp(['Input Correlation = ' num2str(corrIn(2)) char(10) 'Reference Correlation = ' num2str(corrRef(2)) char(10) 'Output Correlation = ' num2str(corrOut(2))])
%     close all
%     scatter(vecRefHist, vecRefHistV, 'linewidth', 2)
%     hold on
%     scatter(vecSim2Bc, vecSim2BcV, [],'red', 'linewidth', 2)
%     scatter(vecSimOut, vecSimOutV, '+')
end     