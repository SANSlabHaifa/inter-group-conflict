function ds = AddCM(ds)
ds.cm.description = 'data transformed to centimeters, with speed and direction of velocity';
ds.cm.t      = cell(ds.raw.numberofsubjects,1);
ds.cm.dx     = cell(ds.raw.numberofsubjects,1);
ds.cm.dy     = cell(ds.raw.numberofsubjects,1);
ds.cm.x      = cell(ds.raw.numberofsubjects,1);
ds.cm.y      = cell(ds.raw.numberofsubjects,1);
ds.cm.dt     = cell(ds.raw.numberofsubjects,1);
ds.cm.v      = cell(ds.raw.numberofsubjects,1);
ds.cm.theta  = cell(ds.raw.numberofsubjects,1);

for subject = 1:ds.raw.numberofsubjects
    ds.cm.t {subject}  =  ds.raw.data                     {subject}(:,1);
    ds.cm.dx{subject}  = (ds.raw.xwidthcm/2) .*ds.raw.data{subject}(:,2);
    ds.cm.x {subject}  = (ds.raw.xwidthcm/2) .*ds.raw.data{subject}(:,4);
    ds.cm.dy{subject}  = (ds.raw.yheightcm/2).*ds.raw.data{subject}(:,3);
    ds.cm.y {subject}  = (ds.raw.yheightcm/2).*ds.raw.data{subject}(:,5);
    
    ds.cm.dt{subject}    = [0; diff(ds.cm.t{subject})];
    ds.cm.v {subject}    = (ds.cm.dx{subject}.^2 + ds.cm.dy{subject}.^2).^(1/2)./ds.cm.dt{subject};
    ds.cm.theta{subject} = atan2(ds.cm.dy{subject},ds.cm.dx{subject});
    
end

return