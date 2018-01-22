function rank = vec_invrank(orgVec, newVec)

rank = nan(size(newVec));

for ii = 1 : numel(newVec)
    rank(ii) = sum(orgVec(:) <= newVec(ii)) / numel(newVec);
end
