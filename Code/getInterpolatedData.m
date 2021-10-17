function ds = getInterpolatedData(ds, old_ver)
ds.int.description = 'time interpolated data';
n = 2;                                       % arbitrary buffer, seconds, for avoiding start and end
startt     = 0;                                % n seconds after start
if old_ver
    endt       = str2num(ds.raw.log{3}(6:end));    % n seconds before herding
else
    endt       = str2num(ds.raw.log{5}(6:end));
end
ds.int.dt  = ds.raw.samplingrateseconds;                  % constant,          used by all subjects
ds.int.t = ((startt+n) :ds.int.dt:(endt-n))';             % same linear range, used by all subjects
ds.int.dx     = cell(ds.raw.numberofsubjects,1);
ds.int.dy     = cell(ds.raw.numberofsubjects,1);
ds.int.x      = cell(ds.raw.numberofsubjects,1);
ds.int.y      = cell(ds.raw.numberofsubjects,1);
ds.int.v      = cell(ds.raw.numberofsubjects,1);
ds.int.theta  = cell(ds.raw.numberofsubjects,1);

lastinterpolationindices = [];
for subject = 1:length(ds.cm.t)
    lastinterpolationindices = [lastinterpolationindices find(ds.cm.t{subject}<endt,1,'last')];
end    
interpolationindices = 1:min(lastinterpolationindices); % limiting range to avoid NaN near end
% disp(lastinterpolationindices)
% disp(interpolationindices)
for j = 1: ds.raw.numberofsubjects
    ds.int.x     {j} = interp1(ds.cm.t{j}(interpolationindices),ds.cm.x    {j}(interpolationindices),ds.int.t);
    ds.int.y     {j} = interp1(ds.cm.t{j}(interpolationindices),ds.cm.y    {j}(interpolationindices),ds.int.t);
    ds.int.dx    {j} = [0;diff(ds.int.x{j})];
    ds.int.dy    {j} = [0;diff(ds.int.y{j})];    
    ds.int.v     {j} = (ds.int.dx{j}.^2 + ds.int.dy{j}.^2).^(1/2)./ds.int.dt;
    ds.int.theta {j} = atan2(ds.int.dy{j},ds.int.dx{j});
end
return