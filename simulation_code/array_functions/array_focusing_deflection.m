function [opas] = array_focusing_deflection(opas)
% OPAS simulation: Beam focusing, multi-beam synthesis, and deflection
% 11/Nov/2025, Yanwei Huang (mrhuangyw1998@gmail.com)
% Parameters

% Two-beam generation using the central and outer subapertures
multi_beam_enabled = isfield(opas, 'multi_beam_enable') ...
    && opas.multi_beam_enable == 1;

if multi_beam_enabled
    if isfield(opas, 'bessel_like_enable') && opas.bessel_like_enable == 1
        error('Bessel-like and multi-beam modes cannot be enabled simultaneously.');
    end

    if opas.helical_wavefront_enable == 1
        error('Helical-wavefront and multi-beam modes cannot be enabled simultaneously.');
    end

    [opas] = array_multi_beam(opas);
    return;
end


if opas.source_focus ~= 0
    opas.time_delays = -(sqrt(opas.bundle_coord(:,1).^2 + opas.bundle_coord(:,2).^2 + opas.source_focus.^2) - opas.source_focus) ./ opas.speed;
    opas.time_delays = opas.time_delays - min(opas.time_delays); % Note: time delay values must be positive
else
    opas.time_delays = zeros(opas.element_num,1);
end

% Time delay for helical wavefront
[opas] = array_helical_wf(opas);

% time delays for beam deflection 
element_pitch_azi = opas.element_pitch; 
element_pitch_ele = opas.element_pitch*cosd(30); 
phase_step_azi = 2*pi*element_pitch_azi*sind(opas.source_deflect_azi)/(opas.speed/opas.source_f0); % phase step angle, azimuthal (axis-x)
phase_step_ele = 2*pi*element_pitch_ele*sind(opas.source_deflect_ele)/(opas.speed/opas.source_f0); % phase step angle, elevation (axis-y)

% Note: always use the top left element (i=1) as the reference when calculating steering phase delay
for i = 1:opas.element_num 
    % Azimuthal time delay calculation
    opas.time_delays(i) = opas.time_delays(i)+((opas.bundle_coord(i,1)-opas.bundle_coord(1,1))/element_pitch_azi)*phase_step_azi/(2*pi)/opas.source_f0; 
    
    % Elevation time delay calculation
    opas.time_delays(i) = opas.time_delays(i)+((opas.bundle_coord(i,2)-opas.bundle_coord(1,2))/element_pitch_ele)*phase_step_ele/(2*pi)/opas.source_f0;
    
end 

end

