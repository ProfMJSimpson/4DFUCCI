function [reddata,yeldata,gredata,arrdata,deddata,cycdata] = radcalcs(reds,yels,gres,arrs,deads,pbin,dayqruniradarr)

% Find all of the cycling red information 
arridxs = ismember(reds(:,1),arrs(:,1)); % Find all arrested cell indices in the red array
cycidxs = ~arridxs;
cycs = reds(cycidxs,:); % List of cycling red cells

if ~isempty(deads)
    centroid = [mean([reds(:,1) ; yels(:,1) ; gres(:,1) ; deads(:,1)]) mean([reds(:,2) ; yels(:,2) ; gres(:,2) ; deads(:,2)]) mean([reds(:,3) ; yels(:,3) ; gres(:,3) ; deads(:,3)])];
else
    centroid = [mean([reds(:,1) ; yels(:,1) ; gres(:,1)]) mean([reds(:,2) ; yels(:,2) ; gres(:,2)]) mean([reds(:,3) ; yels(:,3) ; gres(:,3)])];
end

redrad = sqrt((reds(:,1) - centroid(1)).^2 + (reds(:,2) - centroid(2)).^2 + (reds(:,3) - centroid(3)).^2);
yelrad = sqrt((yels(:,1) - centroid(1)).^2 + (yels(:,2) - centroid(2)).^2 + (yels(:,3) - centroid(3)).^2);
grerad = sqrt((gres(:,1) - centroid(1)).^2 + (gres(:,2) - centroid(2)).^2 + (gres(:,3) - centroid(3)).^2);
arrgreinds = grerad < dayqruniradarr;
grerad(arrgreinds) = [];
if ~isempty(arrs)
    arrrad = sqrt((reds(:,1) - centroid(1)).^2 + (reds(:,2) - centroid(2)).^2 + (reds(:,3) - centroid(3)).^2);
end
if ~isempty(deads)
    dedrad = sqrt((deads(:,1) - centroid(1)).^2 + (deads(:,2) - centroid(2)).^2 + (deads(:,3) - centroid(3)).^2);
end
cycrad = sqrt((cycs(:,1) - centroid(1)).^2 + (cycs(:,2) - centroid(2)).^2 + (cycs(:,3) - centroid(3)).^2); 

[redcount,edgred] = histcounts(redrad,fliplr(pbin));
[yelcount,edgyel] = histcounts(yelrad,fliplr(pbin));
[grecount,edggre] = histcounts(grerad,fliplr(pbin));
if ~isempty(arrs)
    [arrcount,edgarr] = histcounts(arrrad,fliplr(pbin));
end
if ~isempty(deads)
    [dedcount,edgded] = histcounts(dedrad,fliplr(pbin));
end
[cyccount,edgcyc] = histcounts(cycrad,fliplr(pbin));

reddata = redcount;
yeldata = yelcount;
gredata = grecount;
if ~isempty(arrs)
    arrdata = arrcount;
end
if ~isempty(deads)
    deddata = dedcount;
end
cycdata = cyccount;

if isempty(arrs)
    arrdata = zeros(length(pbin)-1,1);
end
if isempty(deads)
    deddata = zeros(length(pbin)-1,1);
end


end