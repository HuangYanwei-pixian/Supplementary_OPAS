%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% MATLAB code for MOPAUS beamforming simulation and visualization
% For manuscript titled 'A two-dimensional miniaturized all-optical 
% phased-array ultrasound source with agile beamforming and frequency tuning'
% 
% Requirement: k-Wave toolbox for the time-domain simulation of acoustic wave fields.
% 
% Created by Yanwei Huang (h.yanwei@wustl.edu) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clc; 
clear; 
close all; 

addpath('array_functions'); % supporting functions of MOPAUS

%% Define parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ultrasound parameters
c0                     = 1480;                  % speed constant [m/s] 
mopaus.source_f0       = 13e6;                  % ultrasound center frequency [Hz]. 61-core tuning range: 5~27
mopaus.source_amp      = 1e6;                   % initial source pressure [Pa]
mopaus.source_cycles   = 1;                     % number of toneburst cycles 
mopaus.speed           = c0;                    % ultrasound speed [m/s] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fiber bundle array parameters
mopaus.d_core = 125e-6; % diameter of each fiber element (include cladding). Unit: m.
mopaus.d_num = 9; % number of fiber cores in diameter. Odd number only! 

% Beamforming parameters
% focusing
mopaus.source_focus    = 1e-3;                  % focal length [m]

% deflection
mopaus.source_deflect_azi  = 0;                 % azimuthal deflection angle (axis-x) [deg]
mopaus.source_deflect_ele  = 0;                 % elevation deflection angle (axis-y) [deg]

% helical wavefront
mopaus.helical_wavefront_enable = 0;      % set if helical wavefront is enabled, 1 for valid
mopaus.tpl_charge = -3;                    % Topological charge

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-wave simulation parameters
% select which processor to run k-Wave code 
%   1: CPU computing 
%   2: GPU computing 
model           = 2;  
parallel.gpu.enableCUDAForwardCompatibility(true)

% Medium parameters
medium.density  = 1000;     % medium density [kg/m^3]
medium.sound_speed = c0;    % ultrasound speed in medium [m/s]

% Computational parameters
% transducer position
translation     = [0, 0, 0]; 
rotation        = [0, 0, 0]; 

% computatonal grid range
grid_size_x     = 2000e-6; % [m] 
grid_size_y     = 2000e-6; % [m] 
grid_size_z     = 4000e-6; % [m] 

% computational parameters
ppw             = 3;        % number of points per wavelength
t_end           = 3e-6;   % total compute time [s]
cfl             = 0.2;      % CFL number
source_x_offset = 20;       % grid points to offset the source (invalid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step digitization
mopaus.time_delay_digitized = 1;   % set if the time delay should be digitized, 1 for valid
mopaus.time_step = 2;      % minimum time step (unit: ns), depending on the system frequency of FPGA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a fiber bundle array
[mopaus] = array_fiber_bundle_coordinates(mopaus, 1); 

% Phase synthesis
[mopaus] = array_focusing_deflection(mopaus); % Beam focusing and deflection

% Time delay digitization
[mopaus] = array_time_digitize(mopaus, 2);

% path to save pictures
save_path = 'sim_results\';

if ~exist(save_path, 'dir')
    mkdir(save_path);
end

if mopaus.helical_wavefront_enable == 1
    file_name = ['_enum' num2str(mopaus.element_num) '_freq' num2str(mopaus.source_f0/1e6) 'MHz_f' num2str(mopaus.source_focus*1e6) 'um_def' num2str(mopaus.source_deflect_azi) '(azi)_' num2str(mopaus.source_deflect_ele) '(ele)deg_TC' num2str(mopaus.tpl_charge) '.fig'];
else
    file_name = ['_enum' num2str(mopaus.element_num) '_freq' num2str(mopaus.source_f0/1e6) 'MHz_f' num2str(mopaus.source_focus*1e6) 'um_def' num2str(mopaus.source_deflect_azi) '(azi)_' num2str(mopaus.source_deflect_ele) '(ele)deg.fig'];
end

% Save time delays
full_file_name0 = [save_path, 'Time delays', file_name];
saveas(gcf, full_file_name0);

% set a depth to show lateral profile
show_lateral_depth = 1e-3;  % [m]

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-wave simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GRID
% calculate the grid spacing based on the PPW and F0
d_vox = c0 / (ppw * mopaus.source_f0);   % [m]

% compute the size of the grid
Nx = roundEven(grid_size_x / d_vox);
Ny = roundEven(grid_size_y / d_vox);
Nz = roundEven(grid_size_z / d_vox);

% create the computational grid
kgrid = kWaveGrid(Nx, d_vox, Ny, d_vox, Nz, d_vox);

% create the time array
kgrid.makeTime(c0, cfl, t_end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate ultrasound field distribution
% create time varying source signals (one for each physical element)
source_sig = mopaus.source_amp .* toneBurst(1/kgrid.dt, mopaus.source_f0, mopaus.source_cycles, 'SignalOffset', round(mopaus.time_delays_show / kgrid.dt));

% create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

% add disc elements
for ind = 1:mopaus.element_num
    % add array element
    karray.addDiscElement([mopaus.bundle_coord(ind,1), mopaus.bundle_coord(ind,2), kgrid.z_vec(1)], mopaus.element_diameter, rotation);
end

% move the array
karray.setArrayPosition(translation, rotation);

% assign binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% plot the off-grid source mask
% voxelPlot(single(source.p_mask));
% title('Off-grid source mask');

% assign source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);

% SENSOR 
% set sensor mask to record central plane
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(:, Ny/2, :) = 1;

% record the pressure
sensor.record = {'p_max'};


%% SIMULATION AND VISUALIZATION
% set input options
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};

if mopaus.helical_wavefront_enable == 0
% run code
switch model
    case 1
        % MATLAB CPU
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
            input_args{:}, ...
            'DataCast', 'single', ...
            'PlotScale', [-1, 1] * mopaus.source_amp);
    case 2
        % MATLAB GPU
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
            input_args{:}, ...
            'DataCast', 'gpuArray-single', ...
            'PlotScale', [-1, 1] * mopaus.source_amp);
end

% reshape data
p_max = reshape(sensor_data.p_max, Nx, Nz);

%% Visualization

show_fwhm_marker = false;
show_peak_marker = true;

z_axis_mm = 1e3 * (kgrid.z_vec + abs(kgrid.z_vec(1)));
x_axis_mm = 1e3 * kgrid.x_vec;

pressure_MPa = p_max / 1e6;
pressure_norm = p_max / max(p_max(:));
pressure_dB = 20 * log10(max(pressure_MPa, realmin));

% pressure field, original pressure

figure;
imagesc(z_axis_mm, x_axis_mm, pressure_MPa);
xlabel('Z-position (Axial) (mm)');
ylabel('X-position (mm)');
axis image;
title('Ultrasound Field (Azi)');
colormap('jet');
cb = colorbar;
title(cb, '(MPa)');

full_file_name1 = [save_path, 'Ultrasound Field (Azi)', file_name];
saveas(gcf, full_file_name1);

% pressure field, normalized

figure;
imagesc(z_axis_mm, x_axis_mm, pressure_norm);
xlabel('Z-position (Axial) (mm)');
ylabel('X-position (mm)');
axis image;
title('Ultrasound Field (Azi, norm)');
colormap('jet');
cb = colorbar;
title(cb, '(a.u.)');

full_file_name6 = [save_path, 'Ultrasound Field (Azi_norm)', file_name];
saveas(gcf, full_file_name6);

% pressure field, dB

figure;
imagesc(z_axis_mm, x_axis_mm, pressure_dB);
xlabel('Z-position (Axial) (mm)');
ylabel('x-position (mm)');
axis image;
title('Ultrasound Field (Azi, decibel)');
colormap('jet');
cb = colorbar;
title(cb, '(dB)');

full_file_name2 = [save_path, 'Ultrasound Field (Azi_decibel)', file_name];
saveas(gcf, full_file_name2);

% axial profile, only for zero deflection angle

if mopaus.source_deflect_azi == 0 && mopaus.source_deflect_ele == 0

    axial_profile = p_max(round(length(kgrid.x_vec) / 2) + 1, :) / 1e6;

    figure;
    plot(z_axis_mm, axial_profile, 'LineWidth', 1.5);
    xlabel('Z-position (Axial) (mm)');
    ylabel('Relative pressure (MPa)');
    title('Sectional View of Ultrasound Field (Axial)');

    if show_fwhm_marker
        hold on;
        add_fwhm_marker_1d(z_axis_mm, axial_profile, 'mm');
    end

    if show_peak_marker
        hold on;
        add_peak_marker_1d(z_axis_mm, axial_profile, 'mm', 'MPa');
    end

    full_file_name3 = [save_path, 'Profile (Axial)', file_name];
    saveas(gcf, full_file_name3);

    figure;
    axial_profile_norm = axial_profile / max(axial_profile);
    plot(z_axis_mm, axial_profile_norm, 'LineWidth', 1.5);
    ylim([0, 1.2]);
    xlabel('Z-position (Axial) (mm)');
    ylabel('Relative pressure (a.u.)');
    title('Sectional View of Ultrasound Field (Axial, norm)');

    if show_fwhm_marker
        hold on;
        add_fwhm_marker_1d(z_axis_mm, axial_profile_norm, 'mm');
    end

    if show_peak_marker
        hold on;
        add_peak_marker_1d(z_axis_mm, axial_profile_norm, 'mm', 'a.u.');
    end

    full_file_name7 = [save_path, 'Profile (Axial_norm)', file_name];
    saveas(gcf, full_file_name7);
end

% lateral profile at the focal plane

if isempty(show_lateral_depth)
    show_lateral_depth = mopaus.source_focus;
end

if show_lateral_depth > 0 && show_lateral_depth <= grid_size_z
    lateral_depth_idx = round(length(kgrid.z_vec) * show_lateral_depth / grid_size_z) + 1;
    lateral_depth_idx = min(max(lateral_depth_idx, 1), length(kgrid.z_vec));
    lateral_profile = p_max(:, lateral_depth_idx) / 1e6;
else
    warning('Depth out of range!\n');
    lateral_profile = [];
end

if ~isempty(lateral_profile)

    figure;
    plot(x_axis_mm, lateral_profile, 'LineWidth', 1.5);
    xlim([-1e3 * grid_size_x / 2, 1e3 * grid_size_x / 2]);
    xlabel('X-position (mm)');
    ylabel('Relative pressure (MPa)');
    title('Sectional View of Ultrasound Field (Azi)');

    if show_fwhm_marker
        hold on;
        add_fwhm_marker_1d(x_axis_mm, lateral_profile, 'mm');
    end

    if show_peak_marker
        hold on;
        add_peak_marker_1d(x_axis_mm, lateral_profile, 'mm', 'MPa');
    end

    full_file_name4 = [save_path, 'Profile (Azi)', file_name];
    saveas(gcf, full_file_name4);

    figure;
    lateral_profile_norm = lateral_profile / max(lateral_profile);
    plot(x_axis_mm, lateral_profile_norm, 'LineWidth', 1.5);
    xlim([-1e3 * grid_size_x / 2, 1e3 * grid_size_x / 2]);
    ylim([0, 1.2]);
    xlabel('X-position (mm)');
    ylabel('Normalized pressure (a.u.)');
    title('Sectional View of Ultrasound Field (Azi, norm)');

    if show_fwhm_marker
        hold on;
        add_fwhm_marker_1d(x_axis_mm, lateral_profile_norm, 'mm');
    end

    if show_peak_marker
        hold on;
        add_peak_marker_1d(x_axis_mm, lateral_profile_norm, 'mm', 'a.u.');
    end

    full_file_name5 = [save_path, 'Profile (Azi_norm)', file_name];
    saveas(gcf, full_file_name5);
end
end

%% simulate and visualize the helical wavefront
if mopaus.helical_wavefront_enable == 1
% X-Y plane visualization for helical wavefront at show_lateral_depth

if isempty(show_lateral_depth)
    show_lateral_depth = mopaus.source_focus;
end

show_xy_window_size = 0.8e-3;   % side length of displayed square region [m]
wavefront_valid_threshold = 0.05;

if show_lateral_depth <= 0 || show_lateral_depth > grid_size_z
    warning('show_lateral_depth is out of range. X-Y helical visualization is skipped.\n');
else
    lateral_depth_idx_xy = get_depth_index_from_source_plane(kgrid.z_vec, show_lateral_depth);
    actual_lateral_depth_m = kgrid.z_vec(lateral_depth_idx_xy) - kgrid.z_vec(1);

    fprintf('\nRunning additional X-Y plane simulation for helical visualization at z = %.3f mm...\n', ...
        1e3 * actual_lateral_depth_m);

    sensor_xy.mask = zeros(Nx, Ny, Nz);
    sensor_xy.mask(:, :, lateral_depth_idx_xy) = 1;

    sensor_xy.record = {'p_max', 'p'};

    switch model
        case 1
            sensor_data_xy = kspaceFirstOrder3D(kgrid, medium, source, sensor_xy, ...
                input_args{:}, ...
                'DataCast', 'single', ...
                'PlotScale', [-1, 1] * mopaus.source_amp);

        case 2
            sensor_data_xy = kspaceFirstOrder3D(kgrid, medium, source, sensor_xy, ...
                input_args{:}, ...
                'DataCast', 'gpuArray-single', ...
                'PlotScale', [-1, 1] * mopaus.source_amp);
    end

    %%
    p_max_xy_vec = gather_if_needed(sensor_data_xy.p_max);
    p_max_xy = reshape(p_max_xy_vec, Nx, Ny);

    pressure_xy_max = max(p_max_xy(:));
    if pressure_xy_max > 0
        pressure_xy_norm = p_max_xy / pressure_xy_max;
    else
        pressure_xy_norm = zeros(size(p_max_xy));
    end

    x_axis_mm_xy = 1e3 * kgrid.x_vec;
    y_axis_mm_xy = 1e3 * kgrid.y_vec;

    x_idx_plot = get_center_window_indices(kgrid.x_vec, show_xy_window_size);
    y_idx_plot = get_center_window_indices(kgrid.y_vec, show_xy_window_size);

    x_axis_mm_xy_plot = x_axis_mm_xy(x_idx_plot);
    y_axis_mm_xy_plot = y_axis_mm_xy(y_idx_plot);

    pressure_xy_norm_plot = pressure_xy_norm(x_idx_plot, y_idx_plot);

    figure;
    imagesc(x_axis_mm_xy_plot, y_axis_mm_xy_plot, pressure_xy_norm_plot.');
    set(gca, 'YDir', 'normal');
    xlabel('X-position (mm)');
    ylabel('Y-position (mm)');
    axis image;
    title(sprintf('X-Y Pressure Field at z = %.3f mm (norm)', ...
        1e3 * actual_lateral_depth_m));
    colormap('jet');
    cb = colorbar;
    title(cb, '(a.u.)');

    full_file_name_xy_pressure = [save_path, 'Ultrasound Field (XY_norm)', file_name];
    saveas(gcf, full_file_name_xy_pressure);

    p_xy_time = gather_if_needed(sensor_data_xy.p);
    p_xy_time = orient_sensor_time_series(p_xy_time, nnz(sensor_xy.mask), length(kgrid.t_array));

    t_array = gather_if_needed(kgrid.t_array);
    t_array = t_array(:).';

    [~, wavefront_peak_idx] = max(abs(p_xy_time), [], 2);
    wavefront_arrival_time_s = t_array(wavefront_peak_idx);

    wavefront_xy_ns = reshape(wavefront_arrival_time_s(:) * 1e9, Nx, Ny);

    wavefront_valid_mask = pressure_xy_norm >= wavefront_valid_threshold;

    if any(wavefront_valid_mask(:))
        wavefront_xy_ns = wavefront_xy_ns - min(wavefront_xy_ns(wavefront_valid_mask));
    else
        warning('No valid pressure region was found for wavefront visualization.');
    end

    wavefront_xy_ns(~wavefront_valid_mask) = NaN;
    wavefront_xy_ns_plot = wavefront_xy_ns(x_idx_plot, y_idx_plot);

    figure;
    imagesc(x_axis_mm_xy_plot, y_axis_mm_xy_plot, wavefront_xy_ns_plot.');
    set(gca, 'YDir', 'normal');
    xlabel('X-position (mm)');
    ylabel('Y-position (mm)');
    axis image;
    title(sprintf('X-Y Wavefront Delay at z = %.3f mm', ...
        1e3 * actual_lateral_depth_m));
    colormap('parula');
    cb = colorbar;
    title(cb, '(ns)');

    full_file_name_xy_wavefront = [save_path, 'Wavefront delay (XY)', file_name];
    saveas(gcf, full_file_name_xy_wavefront);
end
end

%% Compare_profiles_from_figs.m
% Load multiple k-Wave visualization .fig files, extract a horizontal profile
% at mid Y (X-position) and overlay the profiles in one plot.

compare_profile = 0;

if compare_profile == 1
clc; clear; close all;

sim_dir = fullfile(pwd, 'sim_results');

file_list = {
    'Ultrasound Field (Azi)_enum61_freq27MHz_f0um_def0(azi)_0(ele)deg.fig'    
    'Ultrasound Field (Azi)_enum61_freq27MHz_f4000um_def0(azi)_0(ele)deg.fig'
    'Ultrasound Field (Azi)_enum61_freq27MHz_f6000um_def0(azi)_0(ele)deg.fig'
};

y_fraction = 0.5;
normalize_to_peak = false;

fig_out = figure('Color','w');
ax_out = axes('Parent', fig_out);
hold(ax_out, 'on');
grid(ax_out, 'on');

for k = 1:numel(file_list)
    fname = file_list{k};
    fpath = fullfile(sim_dir, fname);

    fprintf('(%d/%d) %s\n', k, numel(file_list), fname);

    if ~isfile(fpath)
        warning('File not found: %s', fpath);
        continue;
    end

    hFig = [];
    try
        hFig = openfig(fpath, 'invisible');

        axAll = findobj(hFig, 'Type', 'axes');
        if isempty(axAll)
            error('No axes found in figure.');
        end

        axKeep = [];
        for ia = 1:numel(axAll)
            tg = get(axAll(ia), 'Tag');
            if (ischar(tg) || isstring(tg)) && (strcmpi(tg,'Colorbar') || strcmpi(tg,'legend'))
                continue;
            end
            axKeep = [axKeep; axAll(ia)]; %#ok<AGROW>
        end
        if isempty(axKeep)
            axKeep = axAll;
        end

        bestObj = [];
        bestN = 0;

        for ia = 1:numel(axKeep)
            hImg  = findobj(axKeep(ia), 'Type', 'image');
            hSurf = findobj(axKeep(ia), 'Type', 'surface');
            cand  = [hImg; hSurf];

            for ic = 1:numel(cand)
                Ctmp = get(cand(ic), 'CData');
                if isnumeric(Ctmp) && ~isempty(Ctmp)
                    n = numel(Ctmp);
                    if n > bestN
                        bestN = n;
                        bestObj = cand(ic);
                    end
                end
            end
        end

        if isempty(bestObj)
            error('No image/surface with numeric CData found.');
        end

        C = get(bestObj, 'CData');
        if ndims(C) == 3
            C = 0.2989*C(:,:,1) + 0.5870*C(:,:,2) + 0.1140*C(:,:,3);
        end
        C = double(C);

        [Ny, Nx] = size(C);

        xData = get(bestObj, 'XData');
        yData = get(bestObj, 'YData');

        xvec = build_axis_vec(xData, Nx);
        yvec = build_axis_vec(yData, Ny);

        y_target = min(yvec) + y_fraction * (max(yvec) - min(yvec));
        [~, iy] = min(abs(yvec - y_target));

        profile = C(iy, :);

        if numel(xvec) >= 2 && xvec(end) < xvec(1)
            xvec = fliplr(xvec);
            profile = fliplr(profile);
        end

        if normalize_to_peak
            m = max(abs(profile));
            if m > 0
                profile = profile / m;
            end
        end

        plot(ax_out, xvec, profile, 'LineWidth', 1.6, 'DisplayName', fname);

    catch ME
        warning('Failed on %s: %s', fname, ME.message);
    end

    if ~isempty(hFig) && isgraphics(hFig)
        close(hFig);
    end
end

xlabel(ax_out, 'Z-position (Axial, mm)');
ylabel(ax_out, 'Field value (from CData)');
title(ax_out, sprintf('Axial distribution of acoustic pressure (y fraction = %.2f)', y_fraction));
legend(ax_out, 'Interpreter','none', 'Location','best');

end

%% helper functions
function vec = build_axis_vec(data, N)
    if isnumeric(data) && isvector(data)
        if numel(data) == 2
            vec = linspace(data(1), data(2), N);
        elseif numel(data) == N
            vec = data(:).';
        else
            vec = 1:N;
        end
    elseif isnumeric(data) && ~isempty(data)
        if size(data,2) == N
            vec = data(1,:);
        elseif size(data,1) == N
            vec = data(:,1).';
        else
            vec = 1:N;
        end
    else
        vec = 1:N;
    end
end

function add_fwhm_marker_1d(x_data, y_data, x_unit)

    x_data = x_data(:);
    y_data = y_data(:);

    [peak_value, peak_idx] = max(y_data);
    half_max = peak_value / 2;

    left_idx = find(y_data(1:peak_idx) <= half_max, 1, 'last');
    right_idx_rel = find(y_data(peak_idx:end) <= half_max, 1, 'first');

    if isempty(left_idx) || isempty(right_idx_rel)
        warning('FWHM could not be determined.');
        return;
    end

    right_idx = peak_idx + right_idx_rel - 1;

    if left_idx >= peak_idx || right_idx <= peak_idx
        warning('FWHM could not be determined.');
        return;
    end

    x_left = interp_halfmax_position( ...
        x_data(left_idx), x_data(left_idx + 1), ...
        y_data(left_idx), y_data(left_idx + 1), ...
        half_max);

    x_right = interp_halfmax_position( ...
        x_data(right_idx - 1), x_data(right_idx), ...
        y_data(right_idx - 1), y_data(right_idx), ...
        half_max);

    fwhm_value = x_right - x_left;

    plot([x_left, x_right], [half_max, half_max], ...
        'k--', 'LineWidth', 1.5);

    plot([x_left, x_left], [0, half_max], ...
        'k:', 'LineWidth', 1.2);

    plot([x_right, x_right], [0, half_max], ...
        'k:', 'LineWidth', 1.2);

    text(mean([x_left, x_right]), half_max, ...
        sprintf('  FWHM = %.3g %s', fwhm_value, x_unit), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center');
end

function add_peak_marker_1d(x_data, y_data, x_unit, y_unit)

    x_data = x_data(:);
    y_data = y_data(:);

    [peak_value, peak_idx] = max(y_data);
    x_peak = x_data(peak_idx);

    plot(x_peak, peak_value, 'ro', ...
        'MarkerSize', 7, ...
        'LineWidth', 1.5);

    plot([x_peak, x_peak], [0, peak_value], ...
        'r:', ...
        'LineWidth', 1.2);

    text(x_peak, peak_value, ...
        sprintf('  max = %.3g %s\n  x = %.3g %s', ...
        peak_value, y_unit, x_peak, x_unit), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'left', ...
        'Color', 'r');
end

function x_half = interp_halfmax_position(x1, x2, y1, y2, half_max)

    if y2 == y1
        x_half = mean([x1, x2]);
    else
        x_half = x1 + (half_max - y1) * (x2 - x1) / (y2 - y1);
    end
end

function idx = get_depth_index_from_source_plane(z_vec, depth_m)

    z_vec = z_vec(:);
    axial_depth_vec = z_vec - z_vec(1);

    [~, idx] = min(abs(axial_depth_vec - depth_m));
    idx = min(max(idx, 1), length(z_vec));
end

function data_out = gather_if_needed(data_in)

    if isa(data_in, 'gpuArray')
        data_out = gather(data_in);
    else
        data_out = data_in;
    end
end

function p_time = orient_sensor_time_series(p_time, num_sensor_points, num_time_points)

    sz = size(p_time);

    if sz(1) == num_sensor_points && sz(2) == num_time_points
        return;
    elseif sz(1) == num_time_points && sz(2) == num_sensor_points
        p_time = p_time.';
    else
        error('Unexpected sensor time-series size. Expected [%d, %d] or [%d, %d].', ...
            num_sensor_points, num_time_points, num_time_points, num_sensor_points);
    end
end

function idx = get_center_window_indices(axis_vec, window_size_m)

    axis_vec = axis_vec(:);

    if isempty(window_size_m) || window_size_m <= 0
        idx = 1:length(axis_vec);
        return;
    end

    center_pos = 0;
    half_width = window_size_m / 2;

    idx = find(axis_vec >= (center_pos - half_width) & axis_vec <= (center_pos + half_width));

    if isempty(idx)
        [~, idx_center] = min(abs(axis_vec - center_pos));
        idx = idx_center;
    end
end
