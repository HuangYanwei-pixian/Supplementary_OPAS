function [opas] = array_helical_wf(opas)
% OPAS simulation: Helical wavefront generation
% 11/Nov/2025, Yanwei Huang (mrhuangyw1998@gmail.com) 
% Parameters

if opas.helical_wavefront_enable == 1 
    opas.bundle_coord_angle = zeros(opas.element_num,1); % Angle of each bundle element under polar coordinate
    for i = 1:opas.element_num 
        % Calculate the angle of each bundle element, unit: rad
        if isnan(opas.bundle_coord(i,2)/opas.bundle_coord(i,1))
            opas.bundle_coord_angle(i) = 0;
        else
            if opas.bundle_coord(i,1) >= 0
                opas.bundle_coord_angle(i) = atan(opas.bundle_coord(i,2)/opas.bundle_coord(i,1));
            else
                opas.bundle_coord_angle(i) = pi+atan(opas.bundle_coord(i,2)/opas.bundle_coord(i,1));
            end
        end

        % Implement time delay
        opas.time_delays(i) = opas.time_delays(i) - (opas.bundle_coord_angle(i)*opas.tpl_charge/(2*pi))/opas.source_f0; 
    end 
end


