clear all;close all;clc;

%Runs the droplet tower cases in cases/BGAL_cases.xlsx and writes the Aspen
%parameter inputs (Ug, mdot_gas, mdot_liquid, eff_CO2/O2/N2) back to each row.
%Profiles are saved to results/<case>_<stage>.mat
%Close BGAL_cases.xlsx in Excel before running.

%USER INPUTS---------------------------------------------------------------
runCases = {}; %case names to run, e.g. {'Water_BL','IL_BL'}; {} = all
runStages = {'abs','des'}; %'abs' and/or 'des'
writeExcel = true; %write outputs back to the case table
%END USER INPUTS ----------------------------------------------------------

root = fileparts(mfilename('fullpath'));
caseFile = fullfile(root,'cases','BGAL_cases.xlsx');
resultsDir = fullfile(root,'results');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end

cases = readtable(caseFile,'Sheet','Cases','TextType','char','VariableNamingRule','preserve');
header = cases.Properties.VariableNames;
col = @(name) colLetter(find(strcmp(header,name)));

%droplet size distribution (same for all cases)
DropletSize = readmatrix(fullfile(root,'DropletDiameterData.csv')); %droplet diameter data (mm)
[GQMOM_model_w, GQMOM_model_x] = dropletQuadrature(DropletSize);

sel = ismember(cases.stage, runStages);
if ~isempty(runCases), sel = sel & ismember(cases.caseName, runCases); end

y_prev = struct('abs', [], 'des', []); %previous solution used as the next initial guess
key_prev = struct('abs', '', 'des', '');

for k = find(sel)'
    c = table2struct(cases(k,:));
    T = c.T_C+273.15;
    if strcmp(c.sorbent,'water')
        p = props_water(T);
    else
        p = props_IL(c.sorbent,T);
    end

    %reuse the previous solution only for the same sorbent and tower height
    key = sprintf('%s_%g', c.sorbent, c.towerHeight);
    if ~strcmp(key, key_prev.(c.stage)), y_prev.(c.stage) = []; end

    if strcmp(c.stage,'abs')
        [out, y_out] = runAbsorber(c, p, GQMOM_model_w, GQMOM_model_x, y_prev.(c.stage));
        fprintf('%-24s abs  Ug %.4f m/s  CO2 eff %6.2f%% (sheet %5.1f)  gas/liquid capture %.4g/%.4g g/s\n', ...
            c.caseName, out.Ug, out.eff_CO2, c.eff_CO2_sheet, out.CO2_captured, out.CO2_capture_rate);
    else
        [out, y_out] = runDesorber(c, p, GQMOM_model_w, GQMOM_model_x, y_prev.(c.stage));
        fprintf('%-24s des  CO2 eff %6.2f%% (sheet %5.1f)\n', c.caseName, out.eff_CO2, c.eff_CO2_sheet);
    end
    y_prev.(c.stage) = y_out;
    key_prev.(c.stage) = key;

    save(fullfile(resultsDir, sprintf('%s_%s.mat', c.caseName, c.stage)), 'out');

    if writeExcel
        r = k+1; %spreadsheet row (header in row 1)
        if strcmp(c.stage,'abs')
            writematrix([out.Ug, out.mdot_gas, out.mdot_liquid], caseFile, 'Sheet', 'Cases', 'Range', sprintf('%s%d', col('Ug'), r));
        else
            writematrix(out.mdot_liquid, caseFile, 'Sheet', 'Cases', 'Range', sprintf('%s%d', col('mdot_liquid'), r));
        end
        writematrix([out.eff_CO2, out.eff_O2, out.eff_N2], caseFile, 'Sheet', 'Cases', 'Range', sprintf('%s%d', col('eff_CO2'), r));
    end
end

function s = colLetter(n)
%spreadsheet column letter for column number n
s = '';
while n > 0
    m = mod(n-1,26);
    s = [char(65+m) s];
    n = floor((n-m-1)/26);
end
end
