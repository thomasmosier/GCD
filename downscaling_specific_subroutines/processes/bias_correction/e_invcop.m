function ec2Bc = e_invcop(dataRef, data2Bc)

ec2Bc = nan(numel(data2Bc(:,1)),1);

copRef = e_cop(dataRef);

for ii = 1 : numel(data2Bc(:,1))
    dist = nan(numel(copRef), 1);
    for ll = 1 : numel(copRef)
        dist(ll) = sqrt(...
             (dataRef(ll,1) - data2Bc(ii,1)).^2 ...
            +(dataRef(ll,2) - data2Bc(ii,2)).^2);
    end

    [~, indUse] = min(dist);
    
    ec2Bc(ii) = copRef(indUse);
end


% [nPts, nDim] = size(data2Bc);
% 
% if nDim ~= 2
%    error('ecopula:dim','The input must be a row vector of ordered pairs, and thus have two columns.');
% end
% 
% %Find MP empirical copula (based on MH)
% ec2Bc = nan(nPts,1);
% for i = 1 : nPts
%     ec2Bc(i) = sum( (dataRef(:,1) <= data2Bc(i,1)).*(dataRef(:,2) <= data2Bc(i,2)) ) / numel(dataRef(:,1));
% end