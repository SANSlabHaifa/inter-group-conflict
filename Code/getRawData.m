function ds = getRawData(basefilename,xwidthcm,yheightcm,samplingrateseconds, old_ver)
%basecd = 'D:\MATLAB\EXPHILA\DATA\';

ds = struct();
ds.raw.filename = basefilename;
ds.raw.description = 'raw data and parameters';
ds.raw.xwidthcm = xwidthcm; 
ds.raw.yheightcm = yheightcm; 
ds.raw.samplingrateseconds = samplingrateseconds;
ds.raw.columnnames = {'time from start of Phase 3(seconds)',...
    'delta x - change from last measurement','delta y - change from last measurement',...
    'x (in range [-1,1])','y (in range [-1,1])'};
ds.raw.data = {};
ds.raw.log = {};

fnlog = [basefilename,'.log'];
logID = fopen(fnlog,'r+');
% logdata = fscanf(logID,'%s');
% disp(fnlog);
tline = fgets(logID);
i = 1;
while ischar(tline)
%     disp(tline)
    tline = fgets(logID);
    ds.raw.log{i} = tline;
    i = i+1;    
end
fclose(logID);
%%%ds.raw.numberofsubjects = (length(ds.raw.log)-6)/2;
if old_ver
    ds.raw.numberofsubjects = str2num(ds.raw.log{5}(10 : end));
else
    ds.raw.numberofsubjects = str2num(ds.raw.log{7}(10 : end));
end
ds.raw.subjects = 0:ds.raw.numberofsubjects-1;

for i = ds.raw.subjects
    fni = [basefilename,'_p',num2str(i),'.log'];
    logi = load(fni);
    ds.raw.data{i+1}    = logi;    
end    
return

