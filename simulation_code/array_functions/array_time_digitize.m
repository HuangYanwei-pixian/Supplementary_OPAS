function [mopaus] = array_time_digitize(mopaus, show_delay)
% Array control system: Time delay digitization
% 11/Nov/2025, Yanwei Huang (mrhuangyw1998@gmail.com) 
% Parameters
%   show_delay: plot the time delay, 1 for valid 

if nargin < 2
    show_delay = 0;
end

multi_beam_enabled = isfield(mopaus, 'multi_beam_enable') && mopaus.multi_beam_enable == 1;

if min(min(mopaus.time_delays)) < 0
    mopaus.time_delays = mopaus.time_delays - min(min(mopaus.time_delays));
end

mopaus.time_delays_show = mopaus.time_delays - min(min(mopaus.time_delays));

% plot time delays if requested
if show_delay == 1
    figure; 
    stem3(mopaus.bundle_coord(:,1),mopaus.bundle_coord(:,2),mopaus.time_delays_show*1e9); 
    zlabel('Time delay (ns)'); 
    hold on 
    theta = 0:pi/50:2*pi; 
    r_total = mopaus.d_num*mopaus.d_core/2; 
    xc_ext = r_total*cos(theta); 
    yc_ext = r_total*sin(theta); 
    plot(xc_ext,yc_ext, 'r'); 
    plot(0,0,'r*'); 
    hold off 
    if multi_beam_enabled
        title_phase_delay = 'Time delay of fiber array elements (two-beam subaperture mode)';
    elseif mopaus.helical_wavefront_enable == 1
        title_phase_delay = ['Time delay of fiber array elements (focus: ' num2str(1000*mopaus.source_focus) 'mm, TC=' num2str(mopaus.tpl_charge) ')'];
    else
        title_phase_delay = ['Time delay of fiber array elements (focus: ' num2str(1000*mopaus.source_focus) 'mm)'];
    end
    title(title_phase_delay); 
end

% Digitalize the time delays
if mopaus.time_delay_digitized == 1  
    fprintf('Time delay has been digitized to %.2f-ns steps.\n', mopaus.time_step);
    for i = 1:mopaus.element_num
        mopaus.time_delays(i) = round(mopaus.time_delays(i)*1e9/mopaus.time_step)*mopaus.time_step/1e9;
    end
end 

if show_delay == 1 && mopaus.time_delay_digitized == 1  
    figure; 
    stem3(mopaus.bundle_coord(:,1),mopaus.bundle_coord(:,2),mopaus.time_delays*1e9); 
    zlabel('Time delay (ns)'); 
    hold on 
    theta = 0:pi/50:2*pi; 
    r_total = mopaus.d_num*mopaus.d_core/2; 
    xc_ext = r_total*cos(theta); 
    yc_ext = r_total*sin(theta); 
    plot(xc_ext,yc_ext, 'r'); 
    plot(0,0,'r*'); 
    hold off 
    if multi_beam_enabled
        title_phase_delay = 'Digitized time delay of fiber array elements (two-beam subaperture mode)';
    elseif mopaus.helical_wavefront_enable == 1
        title_phase_delay = ['Digitized time delay of fiber array elements (focus: ' num2str(1000*mopaus.source_focus) 'mm, TC=' num2str(mopaus.tpl_charge) ')'];
    else
        title_phase_delay = ['Digitized time delay of fiber array elements (focus: ' num2str(1000*mopaus.source_focus) 'mm)'];
    end
    title(title_phase_delay);
end 

% show delay in 2D mode
if show_delay == 2 && mopaus.time_delay_digitized == 1
    if multi_beam_enabled
        title_phase_delay = 'Digitized time delay of fiber array elements (two-beam subaperture mode)';
    else
        title_phase_delay = 'Digitized time delay of fiber array elements (ns)';
    end

    delay_ns = mopaus.time_delays * 1e9;
    delay_min = min(delay_ns);
    delay_max = max(delay_ns);

    figure;
    hold on
    axis equal

    theta = linspace(0, 2*pi, 200);
    r_core = mopaus.d_core / 2;

    for i = 1:mopaus.bundle_num
        xc = mopaus.bundle_coord(i,1) + r_core*cos(theta);
        yc = mopaus.bundle_coord(i,2) + r_core*sin(theta);

        % filled circle with color determined by delay_ns(i)
        patch(xc, yc, delay_ns(i), ...
            'EdgeColor', 'k', ...
            'LineWidth', 1);

        % show element ID at the center
        % text(mopaus.bundle_coord(i,1), mopaus.bundle_coord(i,2), num2str(mopaus.element_id(i)), ...
        %     'HorizontalAlignment', 'center', ...
        %     'VerticalAlignment', 'middle', ...
        %     'FontSize', 10, ...
        %     'FontWeight', 'bold', ...
        %     'Color', 'k');
    end

    % plot external edge
    r_total = mopaus.d_num * mopaus.d_core / 2;
    xc_ext = r_total * cos(theta);
    yc_ext = r_total * sin(theta);
    plot(xc_ext, yc_ext, 'r', 'LineWidth', 1.5);
    plot(0, 0, 'r*');

    title(title_phase_delay);
    xlim([-0.5*mopaus.d_num*mopaus.d_core, 0.5*mopaus.d_num*mopaus.d_core]);
    ylim([-0.5*mopaus.d_num*mopaus.d_core, 0.5*mopaus.d_num*mopaus.d_core]);

    colormap(parula);
    % clim([0, 300]);
    colorbar;

    hold off
end

% reshuffle the order of array elements
mopaus.time_delays_reshuf = zeros(mopaus.bundle_num,1);
for i = 1:mopaus.bundle_num
    mopaus.time_delays_reshuf(i) = mopaus.time_delays(mopaus.order_spi(i));
end
mopaus.time_delays = mopaus.time_delays_reshuf;

mopaus.time_delays_ns = mopaus.time_delays*1e9; 

end

