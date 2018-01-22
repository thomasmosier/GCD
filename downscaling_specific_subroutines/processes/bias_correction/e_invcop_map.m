function map = e_invcop_map(mat1, matRef)

%1st column of 'map' refers to cop1 and 2nd column to cop2 



cop1 = e_invcop(matRef, mat1);
copRef = e_cop(matRef);

    
map = nan(numel(cop1), 2);
for kk = 1 : numel(cop1)
    ind1Curr = find(cop1 == cop1(kk));
    indRefCurr = find(copRef == cop1(kk));
    
    if numel(indRefCurr) == 0
       if  cop1(kk) < min(copRef)
           indRefCurr = find(copRef == min(copRef));
       elseif cop1(kk) > max(copRef)
           indRefCurr = find(copRef == max(copRef));
       else
           error('eInvcopMap:unexpectedValue',['The current empirical copula value is ' ...
               num2str(cop1(kk)) ', which is not contained in the reference copula.']);
       end
    end

    ind1Curr = setdiff(ind1Curr, map(:,1));
    if numel(copRef) >= numel(cop1)
        indRefCurr = setdiff(indRefCurr, map(:,2));
    end

    if numel(ind1Curr) > 1 || numel(indRefCurr) > 1
%             prmsSim = perms(indSimCurr);
%             prmsRef = perms(indRefCurr);

        dist = nan(numel(ind1Curr), numel(indRefCurr));
        for ll = 1 : numel(ind1Curr)
            for mm = 1 : numel(indRefCurr)
                dist(ll,mm) = sqrt(...
                     (mat1(ind1Curr(ll),1) - matRef(indRefCurr(mm),1)).^2 ....
                    +(mat1(ind1Curr(ll),2) - matRef(indRefCurr(mm),2)).^2);
            end
        end
        if all(size(dist)>1)
            [~,ind1Use, indRefUse] = min2d(dist);
        elseif numel(dist(:,1)) == 1
            [~, indRefUse] = min(dist);
            ind1Use = 1;
        else
            [~, ind1Use] = min(dist);
            indRefUse = 1;
        end
        map(kk,:) = [ind1Curr(ind1Use), indRefCurr(indRefUse)]; 
    else
       map(kk,:) = [ind1Curr, indRefCurr]; 
    end
end