% LookNumber - 根据像素间隔计算最佳多视视数
%
% Author: Yinghao Hu
% Repository: https://github.com/Mothanburg/IPSART
function [azimuthLook,rangeLook,resolution] = GetNumLook(azimuthSpacing, rangeSpacing, incidenceAngle, gridSize)

arguments
    azimuthSpacing double
    rangeSpacing double
    incidenceAngle double
    gridSize double {mustBeGreaterThanOrEqual(gridSize, 0)} = 0
end

ground_spacing = rangeSpacing / sind(incidenceAngle);

if azimuthSpacing > ground_spacing
    az_look = 1;
    rg_look = round(azimuthSpacing / ground_spacing);
else
    rg_look = 1;
    az_look = round(ground_spacing / azimuthSpacing);
end

az_res = az_look * azimuthSpacing;
if az_res < gridSize
    az_look = floor(gridSize / azimuthSpacing);
end

rg_res = rg_look * ground_spacing;
if rg_res < gridSize
    rg_look = floor(gridSize / ground_spacing);
end

azimuthLook = az_look;
rangeLook = rg_look;
resolution = [azimuthSpacing * az_look, ground_spacing * rg_look];

end