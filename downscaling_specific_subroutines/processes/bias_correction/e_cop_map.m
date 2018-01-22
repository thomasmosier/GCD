function map = e_cop_map(mat1, mat2)

%1st column of 'map' refers to cop1 and 2nd column to cop2 


cop1 = e_cop(mat1);
cop2 = e_cop(mat2);
    
map = nan(numel(cop1), 2);
for kk = 1 : numel(cop1)
    ind1Curr = find(cop1 == cop1(kk));
    ind2Curr = find(cop2 == cop2(kk));

    ind1Curr = setdiff(ind1Curr, map(:,1));
    ind2Curr = setdiff(ind2Curr, map(:,2));

    if numel(ind1Curr) > 1 || numel(ind2Curr) > 1
%             prmsSim = perms(indSimCurr);
%             prmsRef = perms(indRefCurr);

        dist = nan(numel(ind1Curr), numel(ind2Curr));
        for ll = 1 : numel(ind1Curr)
            for mm = 1 : numel(ind2Curr)
                dist(ll,mm) = sqrt(...
                     (mat1(ind1Curr(ll),1) - mat2(ind2Curr(mm),1)).^2 ....
                    +(mat1(ind1Curr(ll),2) - mat2(ind2Curr(mm),2)).^2);
            end
        end
        if all(size(dist)>1)
            [~,ind1Use, ind2Use] = min2d(dist);
        elseif numel(dist(:,1)) == 1
            [~, ind2Use] = min(dist);
            ind1Use = 1;
        else
            [~, ind1Use] = min(dist);
            ind2Use = 1;
        end
        map(kk,:) = [ind1Curr(ind1Use), ind2Curr(ind2Use)]; 
    else
       map(kk,:) = [ind1Curr, ind2Curr]; 
    end
end