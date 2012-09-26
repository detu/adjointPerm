function well_ind = findEclipseNumberOfWells(well_data,W);
% find the well numbering of the eclipse file of the given wells
wellnames = {W.name}
for i=1:numel(wellnames)
    well_ind(i) = find(strcmp(well_data.WNAME,wellnames{i}));
end