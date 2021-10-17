function ds = AddDistances(ds)
ds.dist.description = 'Center of mass, center of others and distances';

numsubjects = ds.raw.numberofsubjects;

% x and y trajectories of the group's center of mass
xc = zeros(size(ds.int.x{1})); %initailizing 
yc = zeros(size(ds.int.y{1}));
for subject1 = 1:numsubjects
    xc = xc+ds.int.x{subject1};
    yc = yc+ds.int.y{subject1};    
end
ds.dist.centerx = xc./numsubjects;
ds.dist.centery = yc./numsubjects;

% x and y trajectories of the rest of players
ds.dist.othersx = cell(numsubjects,1);
ds.dist.othersy = cell(numsubjects,1);

for subject1 = 1:numsubjects
xo = zeros(size(ds.int.x{1}));
yo = zeros(size(ds.int.y{1}));
for subject2 = [1:subject1-1, subject1+1:numsubjects]
    xo = xo+ds.int.x{subject2};
    yo = yo+ds.int.y{subject2};    
end
ds.dist.othersx{subject1} = xo./(numsubjects-1);
ds.dist.othersy{subject1} = yo./(numsubjects-1);
end


ds.dist.distpair    = cell(numsubjects,numsubjects);
ds.dist.dist2center = cell(numsubjects,1);
ds.dist.dist2others = cell(numsubjects,1);


% distance profiles for each pair
for subject1 = 1:numsubjects
    for subject2 = 1:numsubjects
        x1 = ds.int.x{subject1};
        y1 = ds.int.y{subject1};
        x2 = ds.int.x{subject2};
        y2 = ds.int.y{subject2};
        ds.dist.dist{subject1,subject2} = sqrt((x1-x2).^2 + (y1-y2).^2);
    end
end


% distance profiles from center of mass
for subject1 = 1:numsubjects
    x1 = ds.int.x{subject1};
    y1 = ds.int.y{subject1};
    x2 = ds.dist.centerx;
    y2 = ds.dist.centery;
    ds.dist.dist2center{subject1,1} = sqrt((x1-x2).^2 + (y1-y2).^2);    
end

% distance profiles from center of others
for subject1 = 1:numsubjects
    x1 = ds.int.x{subject1};
    y1 = ds.int.y{subject1};
    x2 = ds.dist.othersx{subject1};
    y2 = ds.dist.othersy{subject1};
    ds.dist.dist2others{subject1,1} = sqrt((x1-x2).^2 + (y1-y2).^2);    
end

% mean distance profiles for each pair
for subject1 = 1:numsubjects
    for subject2 = 1:numsubjects
        ds.dist.meandist{subject1,subject2} = nanmean(ds.dist.dist{subject1,subject2});
        ds.dist.diststd{subject1,subject2} = nanstd(ds.dist.dist{subject1,subject2});
    end
end

% mean distance profiles from center of mass and center of mass with others
for subject1 = 1:numsubjects
    ds.dist.meandist2center{subject1,1} = nanmean(ds.dist.dist2center{subject1,1});    
    ds.dist.meandist2others{subject1,1} = nanmean(ds.dist.dist2others{subject1,1}); 
    ds.dist.dist2centerstd{subject1,1} = nanstd(ds.dist.dist2center{subject1,1});
    sum_dist2others_std = 0;
    for subject2 = [1 : subject1 - 1 subject1 + 1 : numsubjects]
        sum_dist2others_std = sum_dist2others_std + ds.dist.diststd{subject1,subject2};
    end
    ds.dist.dist2othersstd{subject1,1} = sum_dist2others_std / (numsubjects - 1);
end
    
return

