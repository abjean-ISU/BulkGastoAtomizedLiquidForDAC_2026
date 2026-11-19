clear all; close all; clc;

%IL screening figures (same plots as Legacy/ComparisonPlotting.m)
%reads results/screen_<sorbent>_<Ug>mps_abs.mat written by runAll.m

%% Load all results
root = fileparts(mfilename('fullpath'));
names = {'water', 'bmimTf2N', 'emimTf2N', 'bmimDCA', 'bmimSCN'};
Ug_list = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]; %m/s

for k = 1:length(names)
    for j = 1:length(Ug_list)
        out = load(fullfile(root,'results',sprintf('screen_%s_%gmps_abs.mat', names{k}, Ug_list(j)))).out;
        IL(k).Ug_range(j) = out.Ug;
        IL(k).CO2_captured(1,j) = out.CO2_captured; %g/s
        IL(k).capture_eff(1,j) = out.eff_CO2; %%
        IL(k).CO2_purity_abs(1,j) = out.CO2_purity; %%
        IL(k).c_CO2_profile(:,j) = out.c_CO2_x; %g/L, rows = height, cols = Ug
        IL(k).xspan = out.xspan;
    end
end

n_IL = length(IL);

% color scheme - colorblind friendly
colors  = {'#4477AA', '#44AA99', '#CC6677', '#DDAA33', '#AA3377'};
styles  = {'-', '--', '-.', ':', '-'};
markers = {'o', 's', '^', 'd', 'v'};

Ug_range = IL(1).Ug_range;   % [1x6]
xspan    = IL(1).xspan';

IL_names = {'Water', '[bmim][Tf2N]', '[emim][Tf2N]', '[bmim][DCA]', '[bmim][SCN]'};

%% Figure 1 - CO2 Capture Rate
figure('color','white')
hold on
ax=gca;
for k = 1:n_IL
    semilogx(Ug_range, IL(k).CO2_captured(1,:), ...
        [markers{k} styles{k}], ...
        'LineWidth', 2, 'Color', colors{k}, ...
        'DisplayName', IL_names{k})
end
% set log scale AFTER all plotting
set(ax, 'XScale', 'log')
xlabel('Gas Superficial Velocity, U_g (m/s)')
ylabel('CO_2 Capture Rate (g/s)')
legend('Location','northwest')
fontsize(24,'points')

%% Figure 2 - CO2 Capture Efficiency
figure('color','white')
ax=gca;
hold on
for k = 1:n_IL
    semilogx(Ug_range, IL(k).capture_eff(1,:), ...
        [markers{k} styles{k}], ...
        'LineWidth', 2, 'Color', colors{k}, ...
        'DisplayName', IL_names{k})
end
% set log scale AFTER all plotting
set(ax, 'XScale', 'log')
xlabel('Gas Superficial Velocity, U_g (m/s)')
ylabel('CO_2 Capture Efficiency (%)')
legend('Location','northeast')
fontsize(24,'points')

%% Figure 3 - Absorbed CO2 Purity
figure('color','white')
ax = gca;
hold on
for k = 1:n_IL
    semilogx(Ug_range, IL(k).CO2_purity_abs(1,:), ...
        [markers{k} styles{k}], ...
        'LineWidth', 2, 'Color', colors{k}, ...
        'DisplayName', IL_names{k})
end
% set log scale AFTER all plotting
set(ax, 'XScale', 'log')
xlabel('Gas Superficial Velocity, U_g (m/s)')
ylabel('CO_2 Purity in Dissolved Phase (%)')
legend('Location','southeast')
fontsize(24,'points')

%% Index of Ug closest to 0.05 m/s (used by Figure 4 and Figure 5)
[~, j_opt] = min(abs(Ug_range - 0.05));

%% Figure 4 - Dissolved CO2 concentration profiles at Ug=0.05 m/s
figure('color','white')
hold on
for k = 1:n_IL
    profile = IL(k).c_CO2_profile;
    if isvector(profile)
        % already a single profile - just make sure it's a row vector
        yprofile = profile(:).';
    else
        % profile stored per Ug value (rows = height points, cols = Ug_range)
        % pull out the column closest to Ug = 0.05 m/s
        yprofile = profile(:, j_opt).';
    end

    if length(yprofile) ~= length(xspan)
        error('ComparisonPlotting:profileSizeMismatch', ...
            ['c_CO2_profile for %s has %d points but xspan has %d. ', ...
             'Check how c_CO2_profile is stored/sized in results_*.mat.'], ...
            IL_names{k}, length(yprofile), length(xspan))
    end

    plot(xspan, yprofile, ...
        styles{k}, 'LineWidth', 2, 'Color', colors{k}, ...
        'DisplayName', IL_names{k})
end
xlabel('Tower Height (m)')
ylabel('Dissolved CO_2 Concentration (g/L)')
legend('Location','northwest')
fontsize(24,'points')

%% Figure 5 - Bar chart summary at Ug=0.05 m/s (horizontal bars, stacked as rows)
% j_opt (index of Ug closest to 0.05 m/s) was already computed before Figure 4

capture_rates = zeros(1, n_IL);
purities_abs  = zeros(1, n_IL);
efficiencies  = zeros(1, n_IL);

for k = 1:n_IL
    capture_rates(k) = IL(k).CO2_captured(1, j_opt);
    purities_abs(k)  = IL(k).CO2_purity_abs(1, j_opt);
    efficiencies(k)  = IL(k).capture_eff(1, j_opt);
end

rgb = hex2rgb(colors);

% Texture per IL: water stays solid, the ionic liquids each get a distinct
% hatch pattern so bars stay distinguishable in grayscale/print too.
hatch_types = {'none', 'vertical', 'horizontal', 'dots', 'cross'};
hatch_color = [0.15 0.15 0.15];
bar_height  = 0.65;

metric_data   = {capture_rates, efficiencies, purities_abs};
metric_xlabel = {'CO_2 Capture Rate (g/s)', 'Capture Efficiency (%)', 'CO_2 Purity in Dissolved Phase (%)'};
panel_label   = {'(a)', '(b)', '(c)'};

figure('color','white','Position',[100 100 700 900])
for m = 1:3
    subplot(3,1,m)
    hold on
    values = metric_data{m};
    for k = 1:n_IL
        draw_hatched_barh(gca, k, bar_height, values(k), rgb(k,:), hatch_types{k}, hatch_color)
    end
    set(gca, 'YTick', 1:n_IL, 'YTickLabel', IL_names, 'YDir', 'reverse')
    xlabel(metric_xlabel{m})
    box on
    xl = xlim; xlim([0, xl(2)*1.08])   % headroom so bars don't butt against the edge
    text(gca, 0, 1.08, panel_label{m}, ...
        'Units', 'normalized', 'Clipping', 'off', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'FontWeight', 'bold', 'FontSize', 16)
end
fontsize(14,'points')

%% Helper function - hex to rgb
function rgb = hex2rgb(hex_colors)
rgb = zeros(length(hex_colors), 3);
for i = 1:length(hex_colors)
    h = strrep(hex_colors{i}, '#', '');
    rgb(i,:) = [hex2dec(h(1:2)), hex2dec(h(3:4)), hex2dec(h(5:6))] / 255;
end
end

%% Helper function - draw one horizontal bar with an optional hatch texture
function draw_hatched_barh(ax, y_center, height, value, facecolor, hatch_type, hatch_color)
% Draws a single horizontal bar (rectangle patch) on axes ax spanning
% x = [0, value] and y = [y_center-height/2, y_center+height/2], filled
% with facecolor, then overlays a simple hatch pattern on top unless
% hatch_type is 'none'. Patterns are drawn axis-aligned (vertical/
% horizontal lines, dots, cross-hatch) so they look correct regardless
% of the axes' aspect ratio or the bar's length.

y0 = y_center - height/2;
y1 = y_center + height/2;

if value < 0
    x0 = value; x1 = 0;
else
    x0 = 0; x1 = value;
end

% base bar
patch(ax, [x0 x1 x1 x0], [y0 y0 y1 y1], facecolor, ...
    'EdgeColor', 'k', 'LineWidth', 1)

if strcmp(hatch_type, 'none') || (x1 - x0) == 0
    return
end

n_v = 7;   % vertical hatch lines across the bar
n_h = 5;   % horizontal hatch lines across the bar

switch hatch_type
    case 'vertical'
        xs = linspace(x0, x1, n_v + 2); xs = xs(2:end-1);
        for xi = xs
            line(ax, [xi xi], [y0 y1], 'Color', hatch_color, 'LineWidth', 1)
        end
    case 'horizontal'
        ys = linspace(y0, y1, n_h); ys = ys(2:end-1);
        for yi = ys
            line(ax, [x0 x1], [yi yi], 'Color', hatch_color, 'LineWidth', 1)
        end
    case 'dots'
        xs = linspace(x0, x1, n_v + 2); xs = xs(2:end-1);
        ys = linspace(y0, y1, 4); ys = ys(2:end-1);
        [XG, YG] = meshgrid(xs, ys);
        plot(ax, XG(:), YG(:), '.', 'Color', hatch_color, 'MarkerSize', 10)
    case 'cross'
        xs = linspace(x0, x1, n_v + 2); xs = xs(2:end-1);
        for xi = xs
            line(ax, [xi xi], [y0 y1], 'Color', hatch_color, 'LineWidth', 1)
        end
        ys = linspace(y0, y1, n_h); ys = ys(2:end-1);
        for yi = ys
            line(ax, [x0 x1], [yi yi], 'Color', hatch_color, 'LineWidth', 1)
        end
end
end