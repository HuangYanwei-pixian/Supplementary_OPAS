function [mopaus] = array_helical_wf(mopaus)
% Array control system: Helical wavefront generation
% 11/Nov/2025, Yanwei Huang (mrhuangyw1998@gmail.com) 
% Parameters

if mopaus.helical_wavefront_enable == 1 
    mopaus.bundle_coord_angle = zeros(mopaus.element_num,1); % Angle of each bundle element under polar coordinate
    for i = 1:mopaus.element_num 
        % Calculate the angle of each bundle element, unit: rad
        if isnan(mopaus.bundle_coord(i,2)/mopaus.bundle_coord(i,1))
            mopaus.bundle_coord_angle(i) = 0;
        else
            if mopaus.bundle_coord(i,1) >= 0
                mopaus.bundle_coord_angle(i) = atan(mopaus.bundle_coord(i,2)/mopaus.bundle_coord(i,1));
            else
                mopaus.bundle_coord_angle(i) = pi+atan(mopaus.bundle_coord(i,2)/mopaus.bundle_coord(i,1));
            end
        end

        % Implement time delay
        mopaus.time_delays(i) = mopaus.time_delays(i) - (mopaus.bundle_coord_angle(i)*mopaus.tpl_charge/(2*pi))/mopaus.source_f0; 
    end 
end


