function data=get_data_structure_b1(filename,list,frequency)
%get data structure from carsim b1 model
%input: list of indexes from datapoints structure you want to get data from
%output: data structure containing the following fields:
%time: time [s]
%throttle: throttle input
%brake: master clyinder pressure [MPa]
%swa: steering wheel angle [rad]
%wa: wheel angle [rad] found with linear assumption
%vx: longitudinal velocity [m/s]
%vy: lateral velocity [m/s]
%ax: longitudinal acceleration [m/s^2]
%ay: lateral acceleration [m/s^2]
%yr: yaw rate [rad/s]
%dvxdt: derivative of long velocity [m/s^2]
%dvydt: derivative of lat velocity [m/s^2]
%frequency: rate to sample data (Hz)

data=struct('time',[],'throttle',[],'brake',[],'swa',[],'wa',[],'x',[],'y',[],'psi',[],'vx',[],'vy',[],'ax',[],'ay',[],'yr',[],'dvxdt',[],'dvydt',[],'dyrdt',[]);
load(filename,'datapoints')

for i=1:length(list)
    
%time- interpolate data to frequency

if ~exist('frequency','var')
    frequency=mean(diff(datapoints(list(i)).time));
end

timevec=[datapoints(list(i)).time(1):1/frequency:datapoints(list(i)).time(end)]';


data.time=[data.time;timevec];


%inputs
throttle=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,1),timevec);

data.throttle=[data.throttle;throttle];

brake=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,3),timevec);

data.brake=[data.brake;brake];

swa=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,2),timevec);

data.swa=[data.swa;swa];

%states
x=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,4),timevec);

data.x=[data.x;x];

y=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,5),timevec);

data.y=[data.y;y];

psi=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,6),timevec);

data.psi=[data.psi;psi];


vx=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,7),timevec);

data.vx=[data.vx;vx];

vy=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,9),timevec);

data.vy=[data.vy;vy];

ax=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,8),timevec);

data.ax=[data.ax;ax];

ay=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,10),timevec);

data.ay=[data.ay;ay];

yr=interp1(datapoints(list(i)).time,datapoints(list(i)).signals.values(:,16),timevec);

data.yr=[data.yr;yr];

%uncomment these and comment out above if you want numerically
%differentiated with filtered
data.dvxdt=[data.dvxdt;sgolayfilt([diff(vx)./diff(timevec);diff(vx(end-1:end))/diff(timevec(end-1:end))],3,frequency+mod(frequency+1,2))];

data.dvydt=[data.dvydt;sgolayfilt([diff(vy)./diff(timevec);diff(vy(end-1:end))/diff(timevec(end-1:end))],3,frequency+mod(frequency+1,2))];

data.dyrdt=[data.dyrdt;sgolayfilt([diff(yr)./diff(timevec);diff(yr(end-1:end))/diff(timevec(end-1:end))],3,frequency+mod(frequency+1,2))];

end

%use accelerometer data to find accelerations
% data.dvydt=data.ay-data.yr.*data.vx;
% data.dvxdt=data.ax+data.yr.*data.vy;

%use linear relationship for wheel angle
data.wa=60/360*0.3995*data.swa;

%unit conversion from degrees to radians

data.psi=pi/180*data.psi;

data.yr=pi/180*data.yr;

data.dyrdt=pi/180*data.dyrdt;

data.swa=pi/180*data.swa;

data.wa=pi/180*data.wa;

end