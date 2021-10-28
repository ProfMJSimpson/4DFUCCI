function [reddata,yeldata,gredata,arrdata,cycdata] = radcalcs(reds,yels,gres,arrs,pbin,dayqruniradarr)
%--------- RADCALCS ---------%
% radcalcs.m calculates the radial distributions of the different agent
% types (red, yellow, green, arrested, cycling red). This function takes
% Inputs:
%     reds: An Nr x 3 array of red agent X,Y,Z positions
%     yels: An Ny x 3 array of yellow agent X,Y,Z positions
%     gres: An Ng x 3 array of green agent X,Y,Z positions
%     arrs: An Na x 3 array of arrested agent X,Y,Z positions
%     pbin: The vector describing bin edges for the distance to the periphery
%     dayqruniradarr: The arrested radius at Day q, for the ith simulation
%
% Outputs: 
%     reddata: Red agent counts in bin edges described by pbin.
%     yeldata: Yellow agent counts in bin edges described by pbin.
%     gredata: Green agent counts in bin edges described by pbin.
%     arrdata: Arrested red agent counts in bin edges described by pbin.
%     cycdata: Cycling red agent counts in bin edges described by pbin.

% Find all of the cycling red information 
arridxs = ismember(reds(:,1),arrs(:,1)); % Find all arrested cell indices in the red array
cycidxs = ~arridxs; % Cycling red cells are red cells that are not arrested
cycs = reds(cycidxs,:); % List of cycling red cells

% Calculate the location of the centroid
centroid = [mean([reds(:,1) ; yels(:,1) ; gres(:,1)]) mean([reds(:,2) ; yels(:,2) ; gres(:,2)]) mean([reds(:,3) ; yels(:,3) ; gres(:,3)])];

% Calculate the radial coordinate of each agent
redrad = sqrt((reds(:,1) - centroid(1)).^2 + (reds(:,2) - centroid(2)).^2 + (reds(:,3) - centroid(3)).^2);
yelrad = sqrt((yels(:,1) - centroid(1)).^2 + (yels(:,2) - centroid(2)).^2 + (yels(:,3) - centroid(3)).^2);
grerad = sqrt((gres(:,1) - centroid(1)).^2 + (gres(:,2) - centroid(2)).^2 + (gres(:,3) - centroid(3)).^2);
arrgreinds = grerad < dayqruniradarr; % Semantically, we consider any green cell in this region arrested by proxy. There are very few in this region, anyway.
grerad(arrgreinds) = []; % Remove "arrested" green cells
if ~isempty(arrs) % First, check if arrested cell count is 0
    arrrad = sqrt((reds(:,1) - centroid(1)).^2 + (reds(:,2) - centroid(2)).^2 + (reds(:,3) - centroid(3)).^2);
end
cycrad = sqrt((cycs(:,1) - centroid(1)).^2 + (cycs(:,2) - centroid(2)).^2 + (cycs(:,3) - centroid(3)).^2); 

% Fill in the histogram with bin edges "pbin" with agent counts. Flip pbin
% so that it corresponds as a distance from the centroid (radial
% coordinate), as calculated in lines 19-27 
[redcount,edgred] = histcounts(redrad,fliplr(pbin));
[yelcount,edgyel] = histcounts(yelrad,fliplr(pbin));
[grecount,edggre] = histcounts(grerad,fliplr(pbin));
if ~isempty(arrs)
    [arrcount,edgarr] = histcounts(arrrad,fliplr(pbin));
end
[cyccount,edgcyc] = histcounts(cycrad,fliplr(pbin));

% Prepare to export (rename in terms of variables of interest).
reddata = redcount;
yeldata = yelcount;
gredata = grecount;
if ~isempty(arrs)
    arrdata = arrcount;
end
cycdata = cyccount;

if isempty(arrs)
    arrdata = zeros(length(pbin)-1,1);
end


end