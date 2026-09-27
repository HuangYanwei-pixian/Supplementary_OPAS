function [mopaus] = array_fiber_bundle_coordinates(mopaus, show_placement)
% This function generates the Cartesian coordinates of the cores in a fiber bundle.
% by Yanwei Huang (mrhuangyw1998@gmail.com), 24/Mar/2025
% Revisions
%     12Jan2026 / Yanwei: added spiral indexing for RTBF mode. 
% Parameters: 
%   - mopaus.d_num: number of fiber cores in diameter. Odd number only! 
%   - mopaus.d_core: diameter of each fiber core (include cladding). Unit: m.
%   - mopaus.bundle_num: number of fiber cores in the bundle.
%   - mopaus.bundle_coord: Cartesian coordinates of the cores in a fiber bundle.
%   - mopaus.order_spi: Indices of bundle_coord in spiral numbering order
%   - show_placement: show the placement of array elements, 1 for valid

if nargin < 2
    show_placement = 0;
end

%% Compute number of cores in the bundle
if mod(mopaus.d_num,2) ~= 1
    fprintf('Error: The number of fiber cores in diameter must be odd!\n');
    return;
end

d_num_cache = mopaus.d_num;
mopaus.bundle_num = mopaus.d_num;
N_row1 = (mopaus.d_num+1)/2; % number of cores in the first row; aka half of the number of rows 
while d_num_cache > N_row1
    d_num_cache = d_num_cache - 1;
    mopaus.bundle_num = mopaus.bundle_num + 2*d_num_cache;
end
fprintf('Number of fiber cores in the bundle: %d.\n', mopaus.bundle_num);

%% Compute the coordinates of the fiber cores row by row
mopaus.bundle_coord = zeros(mopaus.bundle_num,2); % Cartesian coordinates of the cores, [x,y] in each element
mopaus.bundle_coord_azi = zeros(mopaus.bundle_num,1); % Number of azimuth steps of the cores
mopaus.bundle_coord_ele = zeros(mopaus.bundle_num,1); % Number of elevation steps of the cores

y_shift = mopaus.d_core*cosd(30); % relative position shift in axis-y between adjacent rows ('pitch')

% compute coordinates
cnt_bundle_num = 0;
for i = 1:mopaus.d_num % i: row of the fiber bundle
    
    % compute the number of fiber cores in the i-th row, denoted as N_row
    N_row = i+N_row1-1;
    if N_row > mopaus.d_num
        N_row = mopaus.d_num-(i-N_row1);
    end
    
    % compute the axis-x position of the initial core in the i-th row
    if mod(N_row,2) == 0 % even number of cores in the i-th row
        x_shift_init = -(N_row/2-0.5)*mopaus.d_core;
    else % odd number of cores in the i-th row
        x_shift_init = -((N_row-1)/2)*mopaus.d_core;
    end

    % compute the position of the cores in the i-th row
    for j = 1:N_row % j: the order of the core in the i-th row
        mopaus.bundle_coord(cnt_bundle_num+j,2) = (N_row1-i)*y_shift; % Coordinate in axis-y (ele)
        mopaus.bundle_coord(cnt_bundle_num+j,1) = x_shift_init+(j-1)*mopaus.d_core;% Coordinate in axis-x (azi)
        % coordinates in azimuth and elevation directions, excluding the pitch.
        mopaus.bundle_coord_ele(cnt_bundle_num+j) = (N_row1-i); 
        mopaus.bundle_coord_azi(cnt_bundle_num+j) = (x_shift_init+(j-1)*mopaus.d_core)/(mopaus.d_core/2);
        % Note: In RTBF mode, the bundle's azimuth coordinates are doubled
        % to avoid float number. Therefore, the azimuth pitch size should be halved. 
    end

    cnt_bundle_num = cnt_bundle_num + N_row;
end

% Define array parameters 
mopaus.element_num      = mopaus.bundle_num;       % number of elements
mopaus.element_pitch    = mopaus.d_core;                % pitch [m]
mopaus.element_diameter = 0.84*mopaus.d_core;           % diameter [m]

%% Spiral indexing: center=1, then spiral outward on hex rings
% Convert (x,y) -> axial hex coordinates (u,v) on a flat-top hex lattice:
%   x = d_core*(u + v/2)
%   y = y_shift*v, where y_shift = d_core*cos(30deg)

x = mopaus.bundle_coord(:,1);
y = mopaus.bundle_coord(:,2);

v = round(y / y_shift);
u = round((x - (mopaus.d_core/2).*v) / mopaus.d_core);

% Hex distance (ring index): dist = max(|u|, |v|, |w|), w = -u-v
w = -u - v;
ring = max([abs(u), abs(v), abs(w)], [], 2);
nRings = max(ring);

% Build map from (u,v) -> original index in bundle_coord
key = arrayfun(@(uu,vv) sprintf('%d,%d',uu,vv), u, v, 'UniformOutput', false);
idxMap = containers.Map(key, num2cell(1:mopaus.bundle_num));

% Generate spiral order in axial coords
% Start at (0,0), then ring k starts at (k,0) and walks 6 sides (length k each)
dirs = [-1  1;
        -1  0;
         0 -1;
         1 -1;
         1  0;
         0  1];

spiral_order = zeros(mopaus.bundle_num, 1);
cnt = 1;

% center
if ~isKey(idxMap, "0,0")
    error("Center core (u,v)=(0,0) not found. Check coordinate mapping.");
end
spiral_order(cnt) = idxMap("0,0");

% rings
for k = 1:nRings
    pos = [k, 0]; % start corner on +x direction
    for s = 1:6
        for step = 1:k
            cnt = cnt + 1;
            kstr = sprintf('%d,%d', pos(1), pos(2));
            if ~isKey(idxMap, kstr)
                error("Missing lattice point at (u,v)=(%d,%d).", pos(1), pos(2));
            end
            spiral_order(cnt) = idxMap(kstr);
            pos = pos + dirs(s,:);
        end
    end
end

% Assign label to each original element index
element_id = zeros(mopaus.bundle_num, 1);
element_id(spiral_order) = (1:mopaus.bundle_num).';

% Store results
mopaus.natural_order    = (1:mopaus.bundle_num).';     % label per element (same ordering as bundle_coord)
mopaus.element_id    = element_id;     % label per element (same ordering as bundle_coord)
mopaus.order_spi  = spiral_order;   % indices of bundle_coord in spiral numbering order

%% plot the placement of the fiber cores
if show_placement == 1
title_figure = ['Placement of fiber cores (' num2str(mopaus.bundle_num) ' elements, total diameter=' num2str(mopaus.d_num*mopaus.d_core*1000000) 'um)'];
figure;
% scatter(mopaus.bundle_coord(:,1),mopaus.bundle_coord(:,2),'.');
hold on
axis equal

theta = 0:pi/50:2*pi;
r_core = mopaus.d_core/2;

for i = 1:mopaus.bundle_num
    xc = mopaus.bundle_coord(i,1)+r_core*cos(theta);
    yc = mopaus.bundle_coord(i,2)+r_core*sin(theta);
    plot(xc,yc,'k');

    % show element id
    text(mopaus.bundle_coord(i,1), mopaus.bundle_coord(i,2), num2str(mopaus.element_id(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontSize',10,'FontWeight','bold');
end

% plot external edge
r_total = mopaus.d_num*mopaus.d_core/2;
xc_ext = r_total*cos(theta);
yc_ext = r_total*sin(theta);
plot(xc_ext,yc_ext, 'r');
plot(0,0,'r*');
title(title_figure);
xlim([-0.5*mopaus.d_num*mopaus.d_core, 0.5*mopaus.d_num*mopaus.d_core]);
ylim([-0.5*mopaus.d_num*mopaus.d_core, 0.5*mopaus.d_num*mopaus.d_core]);
hold off
end

end
