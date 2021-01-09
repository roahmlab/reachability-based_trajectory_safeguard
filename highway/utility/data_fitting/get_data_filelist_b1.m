function data=get_data_filelist_b1(list)
%from carsim b1
%input: list of filenames you want to get data from
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

data=struct('time',[],'throttle',[],'brake',[],'swa',[],'wa',[],'vx',[],'vy',[],'ax',[],'jerk',[],'ay',[],'yr',[],'dvxdt',[],'dvydt',[],'dyrdt',[],'gear_no',[]);
for file_index=1:length(list)
load(list{file_index})
%time
data.time=[data.time;datapoints.time];
g=9.80665;

%inputs
data.throttle=[data.throttle;datapoints.signals.values(:,1)];
data.brake=[data.brake;datapoints.signals.values(:,3)];
data.swa=[data.swa;datapoints.signals.values(:,2)];
%states
data.vx=[data.vx;datapoints.signals.values(:,7)];
data.vy=[data.vy;datapoints.signals.values(:,9)];
data.ax=[data.ax;datapoints.signals.values(:,8)];
data.jerk=[data.jerk;get_derivative(datapoints.signals.values(:,8),datapoints.time)];
data.ay=[data.ay;datapoints.signals.values(:,10)];
data.yr=[data.yr;datapoints.signals.values(:,16)];
data.gear_no=[data.gear_no;datapoints.signals.values(:,36)];
%uncomment these and comment out above if you want numerically
%differentiated with filtered
data.dvxdt=[data.dvxdt;smooth([datapoints.time],get_derivative(datapoints.signals.values(:,7),datapoints.time))];

data.dvydt=[data.dvydt;smooth([datapoints.time],get_derivative(datapoints.signals.values(:,9),datapoints.time))];

data.dyrdt=[data.dyrdt;smooth([datapoints.time],get_derivative(datapoints.signals.values(:,16),datapoints.time))];
disp([num2str(file_index),' of ',num2str(length(list))])
end

%use accelerometer data to find accelerations
% data.dvydt=data.ay-data.yr.*data.vx;
% data.dvxdt=data.ax+data.yr.*data.vy;

data.wa=60/360*0.3995*data.swa;


%unit conversion from degrees to radians

data.yr=data.yr;

data.dyrdt=data.dyrdt;

data.swa=pi/180*data.swa;

data.wa=pi/180*data.wa;

end