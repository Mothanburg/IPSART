function [data,dataInfo] = ImportGaoFen3(prefix, dataName)

arguments
    prefix string
    dataName string
end

import matlab.io.xml.dom.Parser

if ~endsWith(prefix, [filesep, "/"])
    prefix = strcat(prefix, filesep);
end
doc = parseFile(Parser, strcat(prefix, filesep, dataName, '.meta.xml'));


% 基本参数
r_num = doc.getElementsByTagName('width').node(1);
a_num = doc.getElementsByTagName('height').node(1);
dataInfo.RangeSize = str2double(r_num.getTextContent);
dataInfo.AzimuthSize = str2double(a_num.getTextContent);

center_time = doc.getElementsByTagName('CenterTime').node(1);
dataInfo.Time = datetime( ...
    center_time.getTextContent, ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSSS', "TimeZone", "UTC" ...
    );

center = doc.getElementsByTagName('imageinfo').node(1).getElementsByTagName('center').node(1).getElementsByTagName('*');
dataInfo.CenterLatitude = str2double(center.node(1));
dataInfo.CenterLongitude = str2double(center.node(2));


% 图像参数
r_space = doc.getElementsByTagName('widthspace').node(1);
a_space = doc.getElementsByTagName('heightspace').node(1);
dataInfo.Image.RangeSpace = str2double(r_space.getTextContent);
dataInfo.Image.AzimuthSpace = str2double(a_space.getTextContent);

inc_near = str2double(doc.getElementsByTagName('incidenceAngleNearRange').node(1).getTextContent);
inc_far = str2double(doc.getElementsByTagName('incidenceAngleFarRange').node(1).getTextContent);
dataInfo.Image.IncidenceAngle = (inc_near + inc_far) / 2;

% 轨道参数
orbit_records = doc.getElementsByTagName('GPSParam');
for idx = 1:orbit_records.getLength
    cur_record = orbit_records.node(idx).getElementsByTagName('*');
    time_stamp = datetime( ...
        cur_record.node(1).getTextContent, ...
        'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSSS', "TimeZone", "UTC" ...
        );

    X = cur_record.node(2).getTextContent;
    Y = cur_record.node(3).getTextContent;
    Z = cur_record.node(4).getTextContent;
    Vx = cur_record.node(5).getTextContent;
    Vy = cur_record.node(6).getTextContent;
    Vz = cur_record.node(7).getTextContent;

    OrbitState(idx).Time = time_stamp;
    OrbitState(idx).PosX = str2double(X);
    OrbitState(idx).PosY = str2double(Y);
    OrbitState(idx).PosZ = str2double(Z);
    OrbitState(idx).VelX = str2double(Vx);
    OrbitState(idx).VelY = str2double(Vy);
    OrbitState(idx).VelZ = str2double(Vz);
end
dataInfo.OrbitRecords = OrbitState;

% 设备和采集参数
f = doc.getElementsByTagName('RadarCenterFrequency').node(1);
dataInfo.Sensing.Frequency = str2double(f) * 1e9;

prf = doc.getElementsByTagName('prf').node(1);
dataInfo.Sensing.PRF = str2double(prf);

rsf = doc.getElementsByTagName('sampleRate').node(1);
dataInfo.Sensing.RangeSampleFrequency = str2double(rsf.getTextContent) * 1e6;

nearRange = doc.getElementsByTagName('imageinfo').node(1).getElementsByTagName('nearRange').node(1);
dataInfo.Sensing.StartRangeTime = str2double(nearRange.getTextContent) / physconst("LightSpeed");

time_imaging = datetime( ...
    doc.getElementsByTagName('imagingTime').node(1).getElementsByTagName('start').node(1).getTextContent, ...
    'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSSSS', "TimeZone", "UTC" ...
    );
[~,m,s] = hms(time_imaging - OrbitState(1).Time);
dataInfo.Sensing.StartAzimuthTime = m * 60 + s;

% 读取数据
polar_mode = doc.getElementsByTagName('polarMode').node(1).getTextContent; % 极化模式
polar_names = ["HH" "HV" "VH" "VV"];
switch polar_mode
    case 'AHV'
        polar = [1 1 1 1];
    otherwise
        error("ImportGaoFen3:invalidPolarMode", "不支持的极化模式：%s", polar_mode);
end

qualify_values = doc.getElementsByTagName('QualifyValue').node(1); % 辐射校正参数
cal_consts = doc.getElementsByTagName('CalibrationConst').node(1);

for idx = 1:4
    if polar(idx)
        raw = readgeoraster(strcat(prefix, strrep(dataName, polar_mode, polar_names(idx)), ".tiff"));
        i = squeeze(double(raw(:,:,1)));
        q = squeeze(double(raw(:,:,2)));
        data.(polar_names(idx)).Data = i + 1i * q;

        qv = qualify_values.getElementsByTagName(polar_names(idx)).node(1).getTextContent;
        data.(polar_names(idx)).QualifyValue = str2double(qv);

        k = cal_consts.getElementsByTagName(polar_names(idx)).node(1).getTextContent;
        data.(polar_names(idx)).CalibrationConst = str2double(k);
    else
        data.(polar_names(idx)).Data = [];
        data.(polar_names(idx)).QualifyValue = nan;
        data.(polar_names(idx)).CalibrationConst = nan;
    end
end

end