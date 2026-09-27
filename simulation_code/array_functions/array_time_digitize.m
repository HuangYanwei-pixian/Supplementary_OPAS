function [opas] = array_time_digitize(opas, show_delay)
% OPAS simulation: Time delay digitization
% 11/Nov/2025, Yanwei Huang (mrhuangyw1998@gmail.com) 
% Parameters
%   show_delay: plot the time delay, 1 for valid 

if nargin < 2
    show_delay = 0;
end

multi_beam_enabled = isfield(opas, 'multi_beam_enable') && opas.multi_beam_enable == 1;

if min(min(opas.time_delays)) < 0
    opas.time_delays = opas.time_delays - min(min(opas.time_delays));
end

opas.time_delays_show = opas.time_delays - min(min(opas.time_delays));

% plot time delays if requested
if show_delay == 1
    figure; 
    stem3(opas.bundle_coord(:,1),opas.bundle_coord(:,2),opas.time_delays_show*1e9); 
    zlabel('Time delay (ns)'); 
    hold on 
    theta = 0:pi/50:2*pi; 
    r_total = opas.d_num*opas.d_core/2; 
    xc_ext = r_total*cos(theta); 
    yc_ext = r_total*sin(theta); 
    plot(xc_ext,yc_ext, 'r'); 
    plot(0,0,'r*'); 
    hold off 
    if multi_beam_enabled
        title_phase_delay = 'Time delay of fiber array elements (two-beam subaperture mode)';
    elseif opas.helical_wavefront_enable == 1
        title_phase_delay = ['Time delay of fiber array elements (focus: ' num2str(1000*opas.source_focus) 'mm, TC=' num2str(opas.tpl_charge) ')'];
    else
        title_phase_delay = ['Time delay of fiber array elements (focus: ' num2str(1000*opas.source_focus) 'mm)'];
    end
    title(title_phase_delay); 
end

% Digitalize the time delays
if opas.time_delay_digitized == 1  
    fprintf('Time delay has been digitized to %.2f-ns steps.\n', opas.time_step);
    for i = 1:opas.element_num
        opas.time_delays(i) = round(opas.time_delays(i)*1e9/opas.time_step)*opas.time_step/1e9;
    end
end 

if show_delay == 1 && opas.time_delay_digitized == 1  
    figure; 
    stem3(opas.bundle_coord(:,1),opas.bundle_coord(:,2),opas.time_delays*1e9); 
    zlabel('Time delay (ns)'); 
    hold on 
    theta = 0:pi/50:2*pi; 
    r_total = opas.d_num*opas.d_core/2; 
    xc_ext = r_total*cos(theta); 
    yc_ext = r_total*sin(theta); 
    plot(xc_ext,yc_ext, 'r'); 
    plot(0,0,'r*'); 
    hold off 
    if multi_beam_enabled
        title_phase_delay = 'Digitized time delay of fiber array elements (two-beam subaperture mode)';
    elseif opas.helical_wavefront_enable == 1
        title_phase_delay = ['Digitized time delay of fiber array elements (focus: ' num2str(1000*opas.source_focus) 'mm, TC=' num2str(opas.tpl_charge) ')'];
    else
        title_phase_delay = ['Digitized time delay of fiber array elements (focus: ' num2str(1000*opas.source_focus) 'mm)'];
    end
    title(title_phase_delay);
end 

% show delay in 2D mode
if show_delay == 2 && opas.time_delay_digitized == 1
    if multi_beam_enabled
        title_phase_delay = 'Digitized time delay of fiber array elements (two-beam subaperture mode)';
    else
        title_phase_delay = 'Digitized time delay of fiber array elements (ns)';
    end

    delay_ns = opas.time_delays * 1e9;
    delay_min = min(delay_ns);
    delay_max = max(delay_ns);

    figure;
    hold on
    axis equal

    theta = linspace(0, 2*pi, 200);
    r_core = opas.d_core / 2;

    for i = 1:opas.bundle_num
        xc = opas.bundle_coord(i,1) + r_core*cos(theta);
        yc = opas.bundle_coord(i,2) + r_core*sin(theta);

        % filled circle with color determined by delay_ns(i)
        patch(xc, yc, delay_ns(i), ...
            'EdgeColor', 'k', ...
            'LineWidth', 1);

        % show element ID at the center
        % text(opas.bundle_coord(i,1), opas.bundle_coord(i,2), num2str(opas.element_id(i)), ...
        %     'HorizontalAlignment', 'center', ...
        %     'VerticalAlignment', 'middle', ...
        %     'FontSize', 10, ...
        %     'FontWeight', 'bold', ...
        %     'Color', 'k');
    end

    % plot external edge
    r_total = opas.d_num * opas.d_core / 2;
    xc_ext = r_total * cos(theta);
    yc_ext = r_total * sin(theta);
    plot(xc_ext, yc_ext, 'r', 'LineWidth', 1.5);
    plot(0, 0, 'r*');

    title(title_phase_delay);
    xlim([-0.5*opas.d_num*opas.d_core, 0.5*opas.d_num*opas.d_core]);
    ylim([-0.5*opas.d_num*opas.d_core, 0.5*opas.d_num*opas.d_core]);

    colormap(parula);
    % clim([0, 300]);
    colorbar;

    hold off
end

% reshuffle the order of array elements
opas.time_delays_reshuf = zeros(opas.bundle_num,1);
for i = 1:opas.bundle_num
    opas.time_delays_reshuf(i) = opas.time_delays(opas.order_spi(i));
end
opas.time_delays = opas.time_delays_reshuf;

opas.time_delays_ns = opas.time_delays*1e9; 

end

