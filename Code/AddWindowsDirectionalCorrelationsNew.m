% AddWindowsDirectionalCorrelationsNew.m
% Adds directional correlation (cosine) of pairs of subjects
% Based on lag:
%    lags are measured in sample #
%    compare subject1 in time t to subject2 in time t+lag
%
% Two types of measurements:
% 1. Normal directional correlation: dc360 - cos of delta_theta
% 2. Opposition insensitive Corrected directional correlation -  dc180
%    treats opposite directions as identical, cosine of 2*delta_theta
%
% zm - identifier of points of no movement, with output as:
% 1 if both are moving
% 0 if both are not moving
% -1 if 1 not moving, 2 moving
% -2 if 1 moving, 2 not moving


function ds = AddWindowsDirectionalCorrelationsNew(ds)
ds.wdcnew.description = 'Moving window directional correlation (new:different lags, 180 degrees, 0 v identification ';
ds.wdcnew.basewindowrange = -5:5;
ds.wdcnew.lagrange      = -10:10;

numsubjects = ds.raw.numberofsubjects;
lenthetas   = length(ds.int.t);
lenlags     = length(ds.wdcnew.lagrange);


ds.wdcnew.dc360   = cell(numsubjects,numsubjects); % directional correlation - cosine( delta_theta)
ds.wdcnew.dc180   = cell(numsubjects,numsubjects); % directional correlation where opposite sides are equal - cosine (2 * delta_theta)
ds.wdcnew.meandc360   = cell(numsubjects,numsubjects); % mean directional correlation - cosine( delta_theta)
ds.wdcnew.meandc180   = cell(numsubjects,numsubjects); % mean directional correlation where opposite sides are equal - cosine (2 * delta_theta)

ds.wdcnew.zm    = cell(numsubjects,numsubjects); % pair testing that both moved
ds.wdcnew.zmall = ones(length(ds.int.t),lenlags); % testing that all pairs moved 


% correlation of thetas
for subject1 = 1:numsubjects
    for subject2 = 1:numsubjects
        ds.wdcnew.dc360{subject1,subject2} = nan(lenthetas,lenlags);
        ds.wdcnew.dc180{subject1,subject2} = nan(lenthetas,lenlags);
        ds.wdcnew.zm{subject1,subject2} = nan(lenthetas,lenlags);
        for lagi = 1:lenlags
            lag = ds.wdcnew.lagrange(lagi); % compare subject1 in time t to subject2 in time t+lag
            if lag<0
                is1 = 1-lag:lenthetas;
                is2 = 1:lenthetas+lag;
            elseif lag==0
                is1 = 1:lenthetas;
                is2 = 1:lenthetas;
            else
                is1 = 1:lenthetas-lag;
                is2 = 1+lag:lenthetas;
            end
            
            theta1 = ds.int.theta{subject1}(is1);
            theta2 = ds.int.theta{subject2}(is2);
            v1     = ds.int.v{subject1}(is1);
            v2     = ds.int.v{subject2}(is2);
            [ds.wdcnew.dc360{subject1,subject2}(is1,lagi),...
                ds.wdcnew.dc180{subject1,subject2}(is1,lagi),...
                ds.wdcnew.zm{subject1,subject2}(is1,lagi)] = ...
                DirectionalCorrelation(theta1,theta2,v1,v2);
        end
        ds.wdcnew.meandc360{subject1,subject2} = nanmean(ds.wdcnew.dc360{subject1,subject2},1);
        ds.wdcnew.meandc180{subject1,subject2} = nanmean(ds.wdcnew.dc180{subject1,subject2},1);
        ds.wdcnew.zmall = ds.wdcnew.zmall .* (ds.wdcnew.zm{subject1,subject2}>0);
    end
end

return


function [dc360,dc180,movements] = DirectionalCorrelation(theta1,theta2,v1,v2)
if ~(length(theta1)==length(theta2) && length(theta2)==length(v1) && length(v1)==length(v2))
    warning('size missmatch')
    return
end
movements = ones(size(theta1));
movements(v1 == 0 & v2 == 0) =  0; % both not moving
movements(v1 == 0 & v2 ~= 0) = -1; % 1 not moving, 2 moving
movements(v1 ~= 0 & v2 == 0) = -2; % 1 moving, 2 not moving

dc360 = cos(theta1-theta2);
dc360(movements==0)  = nan;%1;
dc360(movements==-1) = nan;%0;
dc360(movements==-2) = nan;%0;

dc180 = cos(2.*(theta1-theta2));
dc180(movements==0)  = nan;%1;
dc180(movements==-1) = nan;%0;
dc180(movements==-2) = nan;%0;

return
