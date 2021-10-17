    function ds = initDataStruct(base_path, group, log_file_id, old_ver, session)

xwidthcm = 14;
ywidthcm = 14;
samplingrateseconds = 0.2;

data_path = [base_path, '\', group];
search_string = ['*_', log_file_id, '.log'];
f = dir([data_path, '\', search_string]);
if isempty(f)
    disp('No log file found')
    ds = [];
else
    base_file_name = f(1).name(1 : strfind(f(1).name, '.log') - 1);
    ds = getRawData([data_path, '\', base_file_name],xwidthcm,ywidthcm,samplingrateseconds, old_ver);
    ds = AddCM(ds);
    ds = getInterpolatedData(ds, old_ver);
    ds = AddWindowsDirectionalCorrelationsNew(ds);
    ds = AddDistances(ds);
end