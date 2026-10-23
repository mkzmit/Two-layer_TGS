function [results, fits, fitFigures, sensitivity, tornadoFigure] = TGSAnalysis1()
    % Workflow:
    % - Select one spot from the calibration file and one nominal grating
    % - Subtract polarity-matched baselines; form (POS - NEG)/2
    % - Locate pump arrival and retain the requested time window
    % - Estimate a shared SAW frequency from the selected repeated runs
    % - Fit each repeat independently, using filtered copies only to obtain
    %   thermal starting values. 
    %
    % The optional damped surface acoustic wave (SAW) adds a fitted 
    % lifetime [s]; its frequency [Hz] is fixed to the shared FFT estimate. 
    % Detector amplitudes and the offset are solved by linear least squares
    % at every nonlinear trial.
    %
    % Requires MATLAB with Optimization Toolbox and Signal Processing 
    % Toolbox. Built in MATLAB2025, confirmed to work with MATLAB2026. 
    % All helper functions are local to this file.

    fitFigures = gobjects(0);
    sensitivity = table();
    tornadoFigure = gobjects(0);
    S = validateSettings(defaultSettings());
    traces = prepareTGSData(S);
        if isempty(traces)
            error("TGS:Data", "No complete calibrated runs remain.");
        end
    traces = selectRunsToEvaluate(traces,S.runsToEvaluate);
    initializationTraces = fftThermalInitialization(traces,S);
    [results, fits] = fitTwoLayerTGSRuns(traces,S,initializationTraces);
    disp(results);
    if S.showPlots
        fitFigures = plotTwoLayerTGSResults(traces,fits,S);
    end
    if S.runSensitivity
        [sensitivity, tornadoFigure] = ... 
            startingValueSensitivity(traces,S,results,initializationTraces);
    end
end

%% Settings and Validation
function S = defaultSettings()
    % Centralize editable input, material, and numerical settings

    %% Files and measurement layout
        S.dataDir = "C:\Users\Maken\OneDrive - Massachusetts Institute of Technology\FFUSars\TGS Data rep-20261006T182917Z-1-001\TGS Data rep\MKZ_coatings\WonNb\2026-08-11"; % path to data repository
        S.calFile = "C:\Users\Maken\OneDrive - Massachusetts Institute of Technology\FFUSars\TGS Data rep-20261006T182917Z-1-001\TGS Data rep\MKZ_coatings\WonNb\2026-08-11\WonNb_spot1_postprocessing.txt"; % path to post-processing file with calibration measurements
        S.filePattern = "*.txt";
        S.fileRegex = "(?<nominal>\d+\.\d+)(?:um)?[_-]" + "(?<location>baseline|spot\d+(?:-baseline)?)-" + "(?<polarity>POS|NEG)-(?<run>\d+)\.txt$";
        S.runsToEvaluate = 3; % Fits earliest N runs; inf fits every run
        S.nominalToEvaluate = []; % Grating Spacing to evaluate [um]; empty selects the smallest available
        S.nHeader = 16; % Number of lines until data in files
        S.tCol = 1; % column of time
        S.yCol = 2; % column of signal
        S.tScale = 1; % scale of file time; 1 ~ seconds

    %% Pump arrival and fitting window
        S.tSearch = [0, 100e-9]; % Absolute arrival-search interval [s]
        S.tBaseline = []; % Known pre-pump interval [s]; empty uses first 10% of search
        S.arrivalSmoothPoints = 7; % Used only for arrival detection
        S.arrivalThresholdSigma = 5;
        S.arrivalConsecutivePoints = 3;
        S.nShift = 0; % Extra pre-pump samples retained; does not shift time zero
        S.tFit = [0, inf]; % Pump-relative data interval [s]
        S.fitStartTime = 1e-9; % Exclude t=0; adjust to measured pump/instrument response
        S.nFitPoints = 800; % Evenly spaced acquisition indices after window selection

    %% Acoustic treatment
        S.includeSAW = true; % Fit one damped SAW with fixed FFT frequency
        S.sawVelocityRange = [1.5e3, 3e3]; % Search interval [m/s]
        S.sawBandFraction = 0.10; % Initialization-only removal half-band / frequency
        S.sawMinPeakRatio = 2;
        S.sawSingleTraceMinPeakRatio = 5; % Detection threshold for one run
        S.sawUniformTolerance = 1e-5; % Relative FFT acquisition-grid step tolerance
        S.sawZeroPadFactor = 4; % FFT interpolation; does not improve true resolution
        S.sawMaskTransitionFraction = 0.25; % Taper at initialization-band edges
        S.sawFFTTruncateFraction = 0.9;

    %% Known film and substrate properties
        % Current System ~ 2.2um Tungsten film on Niobium Substrate
        S.L = 2.2e-6; % Film thickness [m]
        S.Cf = 2.56e6; % Film Volumetric heat capacity [J m^-3 K^-1]
        S.Cs = 2.27e6; % Substrate volumetric heat capacity [J m^-3 K^-1]
        S.Ef = 4.11e11; % Film Young's modulus [Pa]
        S.Es = 1.05e11; % Substrate Young's modulus [Pa]
        S.nuf = 0.28; % Film Poisson ratio
        S.nus = 0.40; % Substrate Poisson ratio 
        S.alphathf = 4.42e-6; % Film Thermal expansion coefficient [K^-1]
        S.alphaths = 7.07e-6; % Substrate thermal expansion coefficient [K^-1]

    %% Thermal parameters: [alpha_f,alpha_s,R]
        % Guesses must lie strictly inside positive bounds; fitting uses log10.
        S.p0 = [6.7e-5, 2.7e-5, 1e-8]; % initial guess [m^2/s], [m^2/s], [m^2 K/W]
        S.pLower = [1.5e-7, 1.5e-7, 1e-14]; % lower bound
        S.pUpper = [1e-3, 1e-3, 1e-6]; % upper bound

    %% Numerical reconstruction and optimizer
        S.Nw = 512; % Angular-frequency samples; must be even
        S.wmin = 2 * pi * 1e3; % Smallest positive angular frequency [rad/s]
        S.wmax = 2 * pi * 1e12; % Largest positive angular frequency [rad/s]
        S.Q0 = 1; % Impulse normalization (system dependent)
        S.maxIter = 300; 
        S.maxEval = 3000;
        S.ftol = 1e-8;
        S.xtol = 1e-8;
        S.gtol = 1e-8;
        S.dx = 2e-4; % Relative finite-difference factor applied to log10 coordinates
        S.multistartCount = 15; % tests N guesses for global minimum
        S.multistartSeed = 1; % each seed supplies same guesses each time
        S.display = "off"; % print progress during fit
        S.boundTolerance = 1e-3; % Near-bound tolerance as fraction of log-search span

    %% Output
        S.showPlots = true; % data fitting 
        S.runSensitivity = false; % Toggle starting-value sensitivity refits and tornado plot
        S.sensitivityStartFactors = [0.01, 10]; % One-at-a-time start perturbations
end

function S = validateSettings(S)
    % Check inputs before file loading or model evaluation
    physicalSearchSpace(S);
    positiveFields = ["L", "Cf", "Cs", "Ef", "Es", "Q0", "tScale"];
    for name = positiveFields
        validateattributes(S.(name), {'numeric'}, ...
            {'real', 'scalar', 'finite', 'positive'}, mfilename, char(name));
    end
    for name = ["nuf", "nus"]
        validateattributes(S.(name), {'numeric'}, ...
            {'real', 'scalar', 'finite', '>', -1, '<', 0.5}, mfilename, char(name));
    end
    for name = ["alphathf", "alphaths"]
        validateattributes(S.(name), {'numeric'}, ...
            {'real', 'scalar', 'finite'}, mfilename, char(name));
    end
    validateattributes(S.nHeader, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    validateattributes(S.tCol, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(S.yCol, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(S.showPlots, {'logical'}, {'scalar'});
    validateattributes(S.includeSAW, {'logical'}, {'scalar'});
    if ~isfolder(S.dataDir)
        error("TGS:DataDirectory", ...
            "Data folder does not exist: %s. Update dataDir and calFile in defaultSettings.", S.dataDir);
    end
    validateattributes(S.nFitPoints, {'numeric'}, {'scalar', 'integer', '>=', 10});
    validateattributes(S.fitStartTime, {'numeric'}, {'scalar', 'finite', 'nonnegative'});
    validateattributes(S.nShift, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    validateattributes(S.arrivalSmoothPoints, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(S.arrivalConsecutivePoints, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(S.arrivalThresholdSigma, {'numeric'}, {'scalar', 'finite', 'positive'});
    validateattributes(S.runSensitivity, {'logical'}, {'scalar'});
    validateattributes(S.multistartCount, {'numeric'}, {'scalar', 'integer', '>=', 1});
    validateattributes(S.multistartSeed, {'numeric'}, ...
        {'scalar', 'integer', 'nonnegative', '<=', double(intmax('uint32'))});
    validateattributes(S.dx, {'numeric'}, {'scalar', 'finite', 'positive'});
    validateattributes(S.Nw, {'numeric'}, {'scalar', 'integer', '>=', 4});
    validateattributes(S.wmin, {'numeric'}, {'scalar', 'finite', 'positive'});
    validateattributes(S.wmax, {'numeric'}, {'scalar', 'finite', '>', S.wmin});
    validateattributes(S.boundTolerance, {'numeric'}, {'scalar', 'finite', 'nonnegative', '<', 0.5});
    if mod(S.Nw, 2) ~= 0
        error("TGS:FourierGrid", "Nw must be even.");
    end
    if numel(S.tFit) ~= 2 || isnan(S.tFit(1)) || isnan(S.tFit(2)) || S.tFit(1) >= S.tFit(2)
        error("TGS:FitWindow", "tFit must contain increasing [start,end] times.");
    end
    if numel(S.tSearch) ~= 2 || any(~isfinite(S.tSearch)) || S.tSearch(1) >= S.tSearch(2)
        error("TGS:Arrival", "tSearch must contain increasing finite times.");
    end
    if ~isempty(S.tBaseline) && (numel(S.tBaseline) ~= 2 || ...
            any(~isfinite(S.tBaseline)) || S.tBaseline(1) >= S.tBaseline(2))
        error("TGS:Arrival", "tBaseline must be empty or contain increasing finite times.");
    end
end

%% Data Preparation
function traces = prepareTGSData(S)
    % Build calibrated, baseline-corrected differential traces
    % One structure is returned for each complete POS/NEG pair at the 
    % selected spot and nominal grating. 
    % t and tFull are pump-relative seconds; 
    % Lambda is the calibrated spacing in meters. 
    % y retains the detector signal units. 

    %% Read and validate the selected calibration file
        calListing = dir(S.calFile);
            if isempty(calListing)
                error("TGS:Calibration", "No post-processing calibration files were found.");
            end
            if numel(calListing) ~= 1
                error("TGS:Calibration", ...
                    "Choose exactly one post-processing calibration file; found %d.", numel(calListing));
            end
        calibrationTable = readtable(fullfile(calListing.folder, calListing.name), ...
            "FileType", "text", ...
            "Delimiter", " ", ...
            "MultipleDelimsAsOne", true, ...
            "ReadVariableNames", true, ...
            "VariableNamingRule", "preserve");
        runName = string(calibrationTable.run_name);
        LambdaUm = calibrationSpacing(calibrationTable);
        selectedSpot = calibrationSpot(runName, calListing.name);
        dataPrefix = calibrationDataPrefix(runName, S.fileRegex);
    %% Build the nominal-to-calibrated grating map
        nominalUm = nan(numel(runName), 1);
        for i = 1:numel(runName)
            nominalUm(i) = calibrationNominal(runName(i));
        end
        cal = containers.Map("KeyType", "char", "ValueType", "double");
        nominalValues = unique(nominalUm(isfinite(nominalUm)));
        for i = 1:numel(nominalValues)
            nominal = nominalValues(i);
            candidates = LambdaUm(nominalUm == nominal & ...
                isfinite(LambdaUm) & LambdaUm > 0);
            if isempty(candidates)
                continue
            end
            cal(key(nominal)) = median(candidates) * 1e-6;
        end
    %% Discover and parse measurement files
        listing = dir(fullfile(S.dataDir, S.filePattern));
        listing = listing(startsWith(string({listing.name}), dataPrefix));
        fileTemplate = struct("nominal", [], "spot", [], "polarity", "", "run", [], "baseline", false, "path", "");
        files = repmat(fileTemplate, numel(listing), 1);
        fileCount = 0;
        for i = 1:numel(listing)
            token = regexp(listing(i).name, S.fileRegex, "names", "once");
            if isempty(token)
                continue
            end
            f = fileTemplate;
            f.nominal = str2double(token.nominal);
            f.polarity = string(token.polarity);
            f.run = str2double(token.run);
            f.baseline = contains(token.location, "baseline");
            spotToken = regexp(char(token.location), "spot(?<number>\d+)", "names", "once");
            f.spot = NaN;
            if ~isempty(spotToken)
                f.spot = str2double(spotToken.number);
            end
            f.path = string(fullfile(listing(i).folder, listing(i).name));
            fileCount = fileCount + 1;
            files(fileCount) = f;
        end
        files = files(1:fileCount);
    %% Pair POS and NEG files for the calibration-selected spot
        measurement = find(~[files.baseline] & [files.spot] == selectedSpot);
        baseline = find([files.baseline] & ...
            (isnan([files.spot]) | [files.spot] == selectedSpot));
        if isempty(measurement)
            error("TGS:SpotSelection", ...
                "No measurement files for spot %d were found in the selected directory.", ...
                selectedSpot);
        end
        availableNominal = unique([files(measurement).nominal]);
        nominal = S.nominalToEvaluate;
        if isempty(nominal)
            nominal = min(availableNominal);
        end
        if ~isscalar(nominal) || ~ismember(nominal, availableNominal)
            error("TGS:GratingSelection", ...
                "Choose nominalToEvaluate from: %s um.", ...
                strjoin(string(availableNominal), ", "));
        end
        measurement = measurement([files(measurement).nominal] == nominal);
        pairKey = strings(numel(measurement), 1);
        for i = 1:numel(measurement)
            f = files(measurement(i));
            pairKey(i) = sprintf("%.12g|%d|%d", f.nominal, f.spot, f.run);
        end
        uniquePair = unique(pairKey, "stable");
        traceTemplate = struct("Lambda", [], "t", [], "y", [], ...
            "tFull", [], "yFull", [], "nominal", [], "spot", [], "run", []);
        traces = repmat(traceTemplate, numel(uniquePair), 1);
        traceCount = 0;
    %% Assemble each calibrated differential run
        for i = 1:numel(uniquePair)
            member = measurement(pairKey == uniquePair(i));
            posIndex = member([files(member).polarity] == "POS");
            negIndex = member([files(member).polarity] == "NEG");
            if numel(posIndex) ~= 1 || numel(negIndex) ~= 1
                warning("TGS:SkippedRun", ...
                    "Skipping pair %s: expected exactly one POS and one NEG file.", uniquePair(i));
                continue
            end
            f = files(posIndex);
            mapKey = key(f.nominal);
            if ~isKey(cal, mapKey)
                warning("TGS:SkippedRun", ...
                    "Skipping run %d: no positive calibration for %.4g um.", f.run, f.nominal);
                continue
            end
            try
                [tbp, ybp] = baselineTrace(f.nominal, "POS", files, baseline, S);
                [tbn, ybn] = baselineTrace(f.nominal, "NEG", files, baseline, S);
                [tp, yp] = readTrace(files(posIndex).path, S);
                [tn, yn] = readTrace(files(negIndex).path, S);
                [tp, yp, ybp] = commonGrid(tp, yp, tbp, ybp);
                [tn, yn, ybn] = commonGrid(tn, yn, tbn, ybn);
                [t, yp, yn] = commonGrid(tp, yp - ybp, tn, yn - ybn);
            catch exception
                if startsWith(exception.identifier, "TGS:") || ...
                        startsWith(exception.identifier, "MATLAB:read")
                    warning("TGS:SkippedRun", ...
                        "Skipping run %d at %.4g um: %s", f.run, f.nominal, exception.message);
                    continue
                end
                rethrow(exception);
            end
            % Both polarities have equal weight after their baselines are removed.
                y = 0.5 * (yp - yn);
                t0 = arrivalTime(t, y, S);
                t = t - t0;
            % nShift retains extra samples without changing physical time zero.
                lowerTime = S.tFit(1);
                if S.nShift > 0
                    arrivalIndex = find(t >= 0, 1);
                    lowerTime = min(lowerTime, t(max(1, arrivalIndex - S.nShift)));
                end
            use = t >= lowerTime & t <= S.tFit(2);
            g = traceTemplate;
            g.Lambda = cal(mapKey);
            g.t = t(use);
            g.y = y(use);
            g.tFull = t;
            g.yFull = y;
            g.nominal = f.nominal;
            g.spot = f.spot;
            g.run = f.run;
            traceCount = traceCount + 1;
            traces(traceCount) = g;
        end

    traces = traces(1:traceCount);
end

function [t, y] = baselineTrace(nominal, polarity, files, baseline, S)
    % Return a polarity-matched baseline at the nominal spacing
    % Use an exact baseline when available; otherwise linearly interpolate
    records = baseline([files(baseline).polarity] == polarity);
    nominalList = unique([files(records).nominal]);
    exact = records([files(records).nominal] == nominal);
        if numel(exact) > 1
            error("TGS:Baseline", "Multiple exact %s baselines exist for %.4g um.", polarity, nominal);
        end
        if isscalar(exact)
            [t, y] = readTrace(files(exact).path, S);
            return
        end
    lowerNominal = max(nominalList(nominalList < nominal));
    upperNominal = min(nominalList(nominalList > nominal));
        if isempty(lowerNominal) || isempty(upperNominal)
            error("TGS:Baseline", ...
                "Nominal spacing %.12g is not bracketed by baseline traces.", nominal);
        end
    lower = records([files(records).nominal] == lowerNominal);
    upper = records([files(records).nominal] == upperNominal);
        if numel(lower) ~= 1 || numel(upper) ~= 1
            error("TGS:Baseline", ...
                "Each baseline spacing must contain one %s trace.", polarity);
        end
    [tl, yl] = readTrace(files(lower).path, S);
    [tu, yu] = readTrace(files(upper).path, S);
    [t, yl, yu] = commonGrid(tl, yl, tu, yu);
    weight = (nominal - lowerNominal) / (upperNominal - lowerNominal);
    y = (1 - weight) * yl + weight * yu;
end

function [t, y] = readTrace(path, S)
    %  Read configured columns and retain finite time/signal pairs
    % Time is converted to seconds, sorted, and made unique 
    % Duplicate-time amplitudes are averaged; at least 16 unique samples 
    % must remain.

    try
        A = readmatrix(path, "NumHeaderLines", S.nHeader);
    catch exception
        error("TGS:ReadTrace", "Cannot read %s: %s", path, exception.message);
    end
    if size(A, 2) < max(S.tCol, S.yCol)
        error("TGS:ReadTrace", "Missing time or signal column in %s.", path);
    end
    
    t = A(:, S.tCol) * S.tScale;
    y = A(:, S.yCol);
    finite = isfinite(t) & isfinite(y);
    t = t(finite);
    y = y(finite);
    [t, ~, group] = unique(t);
    y = accumarray(group, y, [], @mean);
    if numel(t) < 16
        error("TGS:ReadTrace", "Fewer than 16 finite unique samples in %s.", path);
    end
end

function [t, first, second] = commonGrid(tFirst, first, tSecond, second)
    % Align two signals to their overlapping acquisition interval.
    % Retain the coarser input time grid and linearly interpolate both 
    % signals onto it.

    if isequal(tFirst, tSecond)
        t = tFirst;
        return
    end
    lower = max(tFirst(1), tSecond(1));
    upper = min(tFirst(end), tSecond(end));
    if median(diff(tFirst)) >= median(diff(tSecond))
        t = tFirst(tFirst >= lower & tFirst <= upper);
    else
        t = tSecond(tSecond >= lower & tSecond <= upper);
    end
    if numel(t) < 16
        error("TGS:CommonGrid", "Traces have fewer than 16 overlapping samples.");
    end
    first = interp1(tFirst, first, t, "linear");
    second = interp1(tSecond, second, t, "linear");
end

function t0 = arrivalTime(t, y, S)
    % Detect a sustained, sign-independent pump onset [s].
    % Estimate offset and noise from the pre-pump interval, threshold a 
    % temporary moving-average signal, then choose the largest nearby 
    % absolute slope.

    search = find(t >= S.tSearch(1) & t <= S.tSearch(2));
        if numel(search) < 16
            error("TGS:Arrival", "Arrival search contains fewer than 16 samples.");
        end
        if isempty(S.tBaseline)
            baseline = search(1:max(5, floor(0.1 * numel(search))));
        else
            baseline = find(t >= S.tBaseline(1) & t <= S.tBaseline(2));
        end
        if numel(baseline) < 5
            error("TGS:Arrival", "At least five known pre-pump samples are required.");
        end
    offset = median(y(baseline));
    % Gaussian-consistent median absolute deviation estimates baseline noise.
        noise = 1.4826 * median(abs(y(baseline) - offset));
    amplitude = max(abs(y(search) - offset));
    noise = max(noise, eps(max(amplitude, realmin)));
    smoothed = movmean(y, S.arrivalSmoothPoints);
    threshold = S.arrivalThresholdSigma * noise;
    above = abs(smoothed(search) - offset) > threshold;
    sustained = conv(double(above), ones(S.arrivalConsecutivePoints, 1), "valid");
    onset = find(sustained >= S.arrivalConsecutivePoints, 1);
    if isempty(onset)
        error("TGS:Arrival", "No sustained pump onset found; check tSearch/tBaseline.");
    end
    if onset <= numel(baseline) && isempty(S.tBaseline)
        error("TGS:Arrival", "Onset overlaps the assumed baseline; set tBaseline explicitly.");
    end
    halfWindow = max(1, ceil(S.arrivalSmoothPoints / 2));
    local = max(1, onset - halfWindow):min(numel(search), onset + halfWindow);
    derivative = gradient(smoothed, t);
    [~, largestSlope] = max(abs(derivative(search(local))));
    index = search(local(largestSlope));
    t0 = t(index);
end

function k = key(value)
    k = sprintf("%.12g", value);
end

function spacingUm = calibrationSpacing(calibrationTable)
    % Read calibrated grating spacing [um].
    % Either grating_spacing_um or grating_value[um] as table header; can 
    % add more as needed
    variableName = string(calibrationTable.Properties.VariableNames);
    supportedName = ["grating_spacing_um", "grating_value[um]"];
    selected = find(ismember(lower(variableName), lower(supportedName)), 1);
        if isempty(selected)
            error("TGS:Calibration", ...
                "Calibration table must contain grating_spacing_um or grating_value[um].");
        end
    spacingUm = double(calibrationTable.(variableName(selected)));
end

function nominalUm = calibrationNominal(runName)
    % Parse the nominal spacing [um] preceding a spot ID
    % Return NaN when a calibration row cannot supply a nominal spacing
    token = regexp(char(runName), "(?<nominal>\d+(?:\.\d+)?)(?:um)?(?=[_-]spot\d+)","names","once");
        if isempty(token)
            nominalUm = NaN;
        else
            nominalUm = str2double(token.nominal);
        end
end

function prefix = calibrationDataPrefix(runName, fileRegex)
    % Identify the filename family from calibration rows.
    % Every run_name must match fileRegex and share the same leading 
    % prefix. Only measurement filenames with that prefix are considered 
    % for fitting.
    prefixByRow = strings(numel(runName), 1);
    for i = 1:numel(runName)
        fileName = runName(i);
        if ~endsWith(fileName, ".txt", "IgnoreCase", true)
            fileName = fileName + ".txt";
        end
        matchStart = regexp(char(fileName), char(fileRegex), "start", "once");
        if isempty(matchStart)
            error("TGS:CalibrationName", ...
                "Calibration run_name '%s' does not match S.fileRegex.", ...
                runName(i));
        end
        prefixByRow(i) = extractBefore(fileName, matchStart);
    end
    prefixList = unique(prefixByRow);
    if numel(prefixList) ~= 1
        error("TGS:CalibrationName", ...
            "Calibration run_name rows identify more than one measurement filename family.");
    end
    prefix = prefixList;
end

function spot = calibrationSpot(runName, calibrationName)
    % Require calibration rows to identify one measurement spot
    % If the calibration filename also contains a spot ID, it must agree

    spotByRow = nan(numel(runName), 1);

    for i = 1:numel(runName)
        token = regexp(char(runName(i)), ...
            "(?:^|[_-])spot(?<spot>\d+)(?=[_-]|$)", ...
            "names", "once");

        if ~isempty(token)
            spotByRow(i) = str2double(token.spot);
        end
    end

    if isempty(spotByRow) || any(~isfinite(spotByRow))
        error("TGS:CalibrationSpot", ...
            "Every calibration run_name must contain a spot identifier such as '_spot1'.");
    end

    spotList = unique(spotByRow);

    if numel(spotList) ~= 1
        error("TGS:CalibrationSpot", ...
            "Choose calibration input for exactly one spot; found spots: %s.", ...
            strjoin(string(spotList), ", "));
    end

    spot = spotList;

    nameToken = regexp(calibrationName, ...
        "(?:^|[_-])spot(?<spot>\d+)(?=[_-]|$)", "names", "once");

    if ~isempty(nameToken) && str2double(nameToken.spot) ~= spot
        error("TGS:CalibrationSpot", ...
            "Calibration filename selects spot %s, but its run_name rows select spot %d.", ...
            nameToken.spot, spot);
    end
end

function traces = selectRunsToEvaluate(traces, runsToEvaluate)
    % Sort complete runs by run number and retain the first N
    % A value of inf selects every available run. 

    validCount = isnumeric(runsToEvaluate) && isreal(runsToEvaluate) && ...
        isscalar(runsToEvaluate) && runsToEvaluate >= 1 && ...
        (isinf(runsToEvaluate) || runsToEvaluate == fix(runsToEvaluate));
        if ~validCount
            error("TGS:RunCount", ...
                "S.runsToEvaluate must be a positive integer or inf.");
        end

    [~, runOrder] = sort([traces.run]);
    traces = traces(runOrder);
    availableRunCount = numel(traces);

    if isfinite(runsToEvaluate)
        if runsToEvaluate > availableRunCount
            error("TGS:RunCount", ...
                "S.runsToEvaluate requests %d runs, but only %d complete " + ...
                "runs are available.", runsToEvaluate, availableRunCount);
        end
        traces = traces(1:runsToEvaluate);
    end
end

function fitTrace = selectFitSamples(trace, S)
    % Apply the fit window, then select acquisition indices
    % Require positive pump-relative time and at least 10 usable samples
    t = trace.t(:);
    y = trace.y(:);
    use = t > 0 & t >= max(S.tFit(1), S.fitStartTime) & t <= S.tFit(2);
        if nnz(use) < 10
            error("TGS:FitWindow", "The positive-time fit window has fewer than 10 samples.");
        end
    t = t(use);
    y = y(use);
    count = min(numel(t), S.nFitPoints);
    indices = unique(round(linspace(1, numel(t), count)));
    fitTrace = trace;
    fitTrace.t = t(indices);
    fitTrace.y = y(indices);
end

%% FFT Initialization 
function initializationTraces = fftThermalInitialization(traces, S)
    % Estimate a common SAW frequency across repeats.
    % Subtract an approximate thermal background for detection, combine 
    % normalized repeat spectra, and remove a narrow band from copies for 
    % thermal pre-fitting. Input repeats must share one grating.

    if isempty(traces) || numel(unique([traces.Lambda])) ~= 1
        error("TGS:FFTTrace", "Supply repeated runs from exactly one grating.");
    end
    
    validateattributes(S.sawFFTTruncateFraction, {'numeric'}, {'scalar', '>', 0, '<=', 1});
    validateattributes(S.sawBandFraction, {'numeric'}, {'scalar', '>', 0, '<', 1});
    validateattributes(S.sawZeroPadFactor, {'numeric'}, {'scalar', 'integer', '>=', 1});
    validateattributes(S.sawUniformTolerance, {'numeric'}, {'scalar', 'finite', 'nonnegative'});
    validateattributes(S.sawMaskTransitionFraction, {'numeric'}, {'scalar', '>=', 0, '<', 1});

    count = numel(traces);
    times = cell(count, 1);
    signals = cell(count, 1);
    residuals = cell(count, 1);
    dt = zeros(count, 1);
    durations = zeros(count, 1);
        for i = 1:count
            [times{i}, signals{i}, dt(i)] = acquisitionGrid(traces(i), S);
            durations(i) = nnz(times{i} > 0 & times{i} <= S.tFit(2)) * dt(i);
            residuals{i} = thermalBackgroundResidual(times{i}, signals{i}, traces(i).Lambda, S);
        end
    searchBand = frequencySearchBand(traces(1).Lambda, dt, S);
    [frequencyHz, peakRatio, detected] = detectFrequency(times, residuals, dt, searchBand, S);
    halfBandwidthHz = NaN;
        if detected
            halfBandwidthHz = max(S.sawBandFraction * frequencyHz, 4 / min(durations));
            if halfBandwidthHz <= 0 || halfBandwidthHz >= frequencyHz || ...
                    frequencyHz + halfBandwidthHz >= min(1 ./ dt) / 2
                error("TGS:SAWBand", "SAW removal band must lie between zero and Nyquist.");
            end
        end
    template = struct("detected", detected, "frequencyHz", frequencyHz, ...
        "peakRatio", peakRatio, "searchBandHz", searchBand, ...
        "halfBandwidthHz", halfBandwidthHz, "sawRMS", 0);
    diagnostics = repmat(template, count, 1);
    initializationTraces = traces;
    for i = 1:count
        saw = zeros(size(signals{i}));
        if detected
            saw = extractBand(times{i}, residuals{i}, dt(i), frequencyHz, halfBandwidthHz, S);
        end
        atFitSamples = interp1(times{i}, saw, traces(i).t(:), "linear");
        if any(~isfinite(atFitSamples))
            error("TGS:FFTTrace", "Initialization samples lie outside the acquisition grid.");
        end
        diagnostics(i).sawRMS = sqrt(mean(saw(times{i} > 0).^2));
        initializationTraces(i).y = traces(i).y(:) - atFitSamples;
        initializationTraces(i).tFull = times{i};
        initializationTraces(i).yFull = signals{i} - saw;
        initializationTraces(i).fftInitialization = diagnostics(i);
    end
end

function [t, y, dt] = acquisitionGrid(trace, S)
    % Validate the full grid for FFT use.
    % Return column vectors and the median time step dt [s]. 
    % Missing or unequal intervals are rejected using sawUniformTolerance

    if isfield(trace, "tFull") && ~isempty(trace.tFull)
        t = trace.tFull(:);
        y = trace.yFull(:);
    else
        t = trace.t(:);
        y = trace.y(:);
    end
    
    if numel(t) ~= numel(y) || numel(t) < 16 || ...
            any(~isfinite(t)) || any(~isfinite(y)) || any(diff(t) <= 0)
        error("TGS:FFTGrid", "FFT input needs finite, paired, increasing acquisition samples.");
    end
    
    dt = median(diff(t));
    
    if max(abs(diff(t) - dt)) / dt > S.sawUniformTolerance
        error("TGS:FFTGrid", ...
            "FFT acquisition grid is nonuniform. Restore missing samples on a uniform grid before FFT use.");
    end
    
    if nnz(t > 0 & t <= S.tFit(2)) < 16
        error("TGS:FFTGrid", "At least 16 positive-time FFT samples are required.");
    end
end

function residual = thermalBackgroundResidual(t, y, Lambda, S)
    % Remove a coarse physical background for the FFT
    % Evaluate the thermal model at p0 and profile temperature, 
    % displacement, and offset amplitudes. 
    % This residual is used only for initialization.
    use = t > 0 & t <= S.tFit(2);
    [temperature, displacement] = TwoLayerModel(t(use), Lambda, S.p0, S);
    design = [temperature, displacement, ones(nnz(use), 1)];
    background = zeros(size(y));
    background(use) = fitLinearComponents(design, y(use));
        if any(t < 0)
            background(~use) = median(y(t < 0));
        end
    residual = y - background;
end

function band = frequencySearchBand(Lambda, dt, S)
    % Convert velocity bounds [m/s] to frequency bounds [Hz]
    band = S.sawVelocityRange(:).' / Lambda;
        if numel(band) ~= 2 || any(~isfinite(band)) || band(1) <= 0 || band(2) < band(1)
            error("TGS:SAWSearch", "sawVelocityRange must contain positive increasing velocities.");
        end
    band(2) = min(band(2), 0.99 * min(1 ./ dt) / 2);
        if band(1) >= band(2)
            error("TGS:SAWSearch", "SAW search range lies outside the acquisition Nyquist limit.");
        end
end

function [frequencyHz, peakRatio, detected] = detectFrequency(times, residuals, dt, band, S)
    % Locate a shared acoustic peak in normalized repeat spectra
    % Differentiate detrended residuals, apply a Hann taper, and take the 
    % median spectrum across runs. 
    % Refine the highest interior peak by a three-point parabola in log 
    % amplitude 
    % Zero padding interpolates bins 
    % Record duration sets the underlying frequency resolution 
    % No peak returns NaN and false.
    count = numel(times);
    frequencies = cell(count, 1);
    amplitudes = cell(count, 1);
    binSteps = zeros(count, 1);
    for i = 1:count
        postPump = find(times{i} > 0 & times{i} <= S.tFit(2));
        n = min(numel(postPump), max(16, floor(S.sawFFTTruncateFraction * numel(postPump))));
        signal = endpointDetrend(residuals{i}(postPump(1:n)));
        derivative = diff(signal) / dt(i);
        derivative = derivative - mean(derivative);
        taper = 0.5 - 0.5 * cos(2 * pi * (0:numel(derivative) - 1).' / (numel(derivative) - 1));
        nfft = 2^nextpow2(S.sawZeroPadFactor * numel(derivative));
        spectrum = abs(fft(derivative .* taper, nfft));
        f = (0:nfft / 2).' / (nfft * dt(i));
        spectrum = spectrum(1:nfft / 2 + 1);
        use = f >= band(1) & f <= band(2);
        if nnz(use) < 5 || diff(band) < 4 / (n * dt(i))
            error("TGS:SAWSearch", "SAW search band has too few resolvable frequencies.");
        end
        scale = max(median(spectrum(use)), realmin);
        frequencies{i} = f;
        amplitudes{i} = spectrum / scale;
        binSteps(i) = 1 / (nfft * dt(i));
    end
    common = (band(1):min(binSteps):band(2)).';
    amplitude = zeros(numel(common), count);
    for i = 1:count
        amplitude(:, i) = interp1(frequencies{i}, amplitudes{i}, common, "linear");
    end
    aggregate = median(amplitude, 2);
    localPeak = find(aggregate(2:end - 1) > aggregate(1:end - 2) & ...
        aggregate(2:end - 1) >= aggregate(3:end)) + 1;
    frequencyHz = NaN;
    peakRatio = 0;
    detected = false;
    if isempty(localPeak)
        return
    end
    [height, selected] = max(aggregate(localPeak));
    peakRatio = height / max(median(aggregate), realmin);
    threshold = S.sawMinPeakRatio;
    if count == 1
        threshold = max(threshold, S.sawSingleTraceMinPeakRatio);
    end
    detected = peakRatio >= threshold;
    if detected
        peak = localPeak(selected);
        logAmplitude = log(max(aggregate(peak - 1:peak + 1), realmin));
        curvature = logAmplitude(1) - 2 * logAmplitude(2) + logAmplitude(3);
        offset = 0;
        if curvature < 0
            offset = 0.5 * (logAmplitude(1) - logAmplitude(3)) / curvature;
            offset = min(max(offset, -0.5), 0.5);
        end
        frequencyHz = common(peak) + offset * min(binSteps);
    end
end

function saw = extractBand(t, residual, dt, frequencyHz, halfBandwidthHz, S)
    % Reconstruct the SAW band for an initialization-only subtraction.
    postPump = find(t > 0 & t <= S.tFit(2));
    signal = endpointDetrend(residual(postPump));
    n = numel(signal);
    nfft = 2^nextpow2(3 * n);
    padded = zeros(nfft, 1);
    indices = n + (0:n - 1);
    padded(indices) = signal;
    f = (0:nfft - 1).' / (nfft * dt);
    f(f > 1 / (2 * dt)) = f(f > 1 / (2 * dt)) - 1 / dt;
    distance = abs(abs(f) - frequencyHz);
    transitionWidth = S.sawMaskTransitionFraction * halfBandwidthHz;
    coreWidth = halfBandwidthHz - transitionWidth;
    mask = double(distance <= coreWidth);
    if transitionWidth > 0
        edge = distance > coreWidth & distance < halfBandwidthHz;
        mask(edge) = 0.5 * (1 + cos(pi * (distance(edge) - coreWidth) / transitionWidth));
    end
    filtered = real(ifft(fft(padded) .* mask));
    saw = zeros(size(residual));
    saw(postPump) = filtered(indices);
end

function signal = endpointDetrend(signal)
    % Subtract the line joining the first and last sample
    signal = signal(:);
    signal = signal - linspace(signal(1), signal(end), numel(signal)).';
end

%% Nonlinear fitting and linear amplitudes
function [results, fits] = fitTwoLayerTGSRuns(traces, S, initializationTraces)
    % Fit repeats independently at one spot and one grating
    % Match initialization copies to run IDs, retain one fit structure 
    % per run, and assemble physical parameters and diagnostics into the 
    % results table. 
    if isempty(traces) || numel(unique([traces.Lambda])) ~= 1 || ...
            numel(unique([traces.spot])) ~= 1
        error("TGS:RunSelection", "Supply runs from exactly one spot and grating.");
    end
    [~, order] = sort([traces.run]);
    traces = traces(order);
    if numel(initializationTraces) ~= numel(traces)
        error("TGS:Initialization", "Initialization must contain one copy per raw run.");
    end
    [found, initializationOrder] = ismember([traces.run], [initializationTraces.run]);
    if ~all(found) || numel(unique(initializationOrder)) ~= numel(traces)
        error("TGS:Initialization", "Initialization run identifiers do not match raw runs.");
    end
    initializationTraces = initializationTraces(initializationOrder);
    fitCell = cell(numel(traces), 1);
    for i = 1:numel(traces)
        fitCell{i} = fitTwoLayerTGS(traces(i), S, initializationTraces(i));
    end
    fits = vertcat(fitCell{:});
    run = [traces.run].';
    p = vertcat(fits.p);
    bound = vertcat(fits.atBound);
    results = table(run, p(:, 1), p(:, 2), p(:, 3), [fits.selectedStart].', ...
        [fits.converged].', [fits.exitflag].', [fits.rmse].', ...
        bound(:, 1), bound(:, 2), bound(:, 3), [fits.sawActive].', ...
        [fits.sawFrequencyHz].', [fits.sawTau].', [fits.sawAtBound].', ...
        [fits.designCondition].', [fits.objective].', [fits.iterations].', ...
        [fits.functionCount].', ...
        'VariableNames', {'Run', 'AlphaF', 'AlphaS', 'R', 'SelectedStart', ...
        'Converged', 'ExitFlag', 'RMSE', 'AlphaFAtBound', 'AlphaSAtBound', 'RAtBound', ...
        'SAWActive', 'SAWFrequencyHz', 'SAWTau', 'SAWTauAtBound', 'DesignCondition', ...
        'NormalizedResnorm', 'Iterations', 'FunctionCount'});
end

function fit = fitTwoLayerTGS(trace, S, initializationTrace)
    % Fit one prepared, unfiltered run with profiled amplitudes
    % Nonlinear coordinates are log10([alpha_f,alpha_s,R]) 
    % SAW frequency is fixed from the shared FFT estimate
    if ~isscalar(trace)
        error("TGS:OneRunFit", "Use fitTwoLayerTGSRuns for repeated runs.");
    end
    [x0, lower, upper] = physicalSearchSpace(S);
    optimizerOptions = optimoptions("lsqnonlin", ...
        "Algorithm", "trust-region-reflective", ...
        "FiniteDifferenceType", "forward", ...
        "FiniteDifferenceStepSize", S.dx, ...
        "MaxIterations", S.maxIter, "MaxFunctionEvaluations", S.maxEval, ...
        "FunctionTolerance", S.ftol, "StepTolerance", S.xtol, ...
        "OptimalityTolerance", S.gtol, "Display", S.display);
    fitTrace = selectFitSamples(trace, S);
    prefitTrace = selectFitSamples(initializationTrace, S);
    scale = residualScale(fitTrace.y);
    prefitScale = residualScale(prefitTrace.y);
    noSAW = struct("active", false, "frequencyHz", NaN);
    prefitResidual = @(x) (profilePrediction(x, prefitTrace, S, noSAW) - prefitTrace.y) / prefitScale;
    try
        xThermal = lsqnonlin(prefitResidual, x0, lower, upper, optimizerOptions);
    catch exception
        if strcmp(exception.identifier, "MATLAB:OperationTerminated")
            rethrow(exception);
        end
        warning("TGS:Prefit", ...
            "Thermal initialization failed for run %d; using p0. %s", trace.run, exception.message);
        xThermal = x0;
    end
    diagnostic = initializationTrace.fftInitialization;
    saw = struct("active", logical(S.includeSAW) && diagnostic.detected, ...
        "frequencyHz", diagnostic.frequencyHz);
    primaryStarts = [xThermal(:).'; x0];
    if saw.active
        [tau0, tauBounds] = acousticLifetime(fitTrace, saw.frequencyHz);
        primaryStarts(:, 4) = log10(tau0);
        lower(4) = log10(tauBounds(1));
        upper(4) = log10(tauBounds(2));
        if max(diff(fitTrace.t)) * saw.frequencyHz > 0.25
            warning("TGS:SAWSampling", ...
                "Selected samples provide fewer than four points per SAW period. Increase nFitPoints.");
        end
    end
    starts = multistartPoints(primaryStarts, lower, upper, S);
    residualFcn = @(x) (profilePrediction(x, fitTrace, S, saw) - fitTrace.y) / scale;
    trialTemplate = struct("x", nan(1, numel(lower)), "resnorm", Inf, ...
        "exitflag", NaN, "iterations", NaN, "functionCount", NaN, "errorMessage", "");
    trials = repmat(trialTemplate, size(starts, 1), 1);
    selected = NaN;
    bestConverged = false;
    bestResnorm = Inf;
    for i = 1:size(starts, 1)
        try
            [x, resnorm, ~, exitflag, output] = lsqnonlin( residualFcn, starts(i, :), lower, upper, optimizerOptions);
        catch exception
            if strcmp(exception.identifier, "MATLAB:OperationTerminated")
                rethrow(exception);
            end
            trials(i).errorMessage = string(exception.message);
            continue
        end
        trials(i).x = x(:).';
        trials(i).resnorm = resnorm;
        trials(i).exitflag = exitflag;
        trials(i).iterations = output.iterations;
        trials(i).functionCount = output.funcCount;
        finite = isfinite(resnorm) && all(isfinite(x));
        converged = finite && exitflag > 0;
        better = finite && (isnan(selected) || (converged && ~bestConverged) || (converged == bestConverged && resnorm < bestResnorm));
        if better
            selected = i;
            bestResnorm = resnorm;
            bestConverged = converged;
        end
    end
    if isnan(selected)
        messages = string({trials.errorMessage});
        first = find(strlength(messages) > 0, 1);
        detail = "No finite objective was returned.";
        if ~isempty(first)
            detail = messages(first);
        end
        error("TGS:Optimization", "All starts failed. %s", detail);
    end
    if ~bestConverged
        warning("TGS:Optimization", "Selected finite fit did not report convergence.");
    end

    bestX = trials(selected).x;
    [prediction, coefficients, conditionNumber, design] = ...
        profilePrediction(bestX, fitTrace, S, saw);
    nearBound = min(bestX - lower, upper - bestX) <= S.boundTolerance * (upper - lower);
    fit.p = 10.^bestX(1:3);
    fit.run = trace.run;
    fit.Lambda = trace.Lambda;
    fit.selectedStart = selected;
    fit.objective = bestResnorm;
    fit.residualScale = scale;
    fit.converged = bestConverged;
    fit.exitflag = trials(selected).exitflag;
    fit.iterations = trials(selected).iterations;
    fit.functionCount = trials(selected).functionCount;
    fit.atBound = nearBound(1:3);
    fit.t = fitTrace.t;
    fit.y = fitTrace.y;
    fit.yFit = prediction;
    fit.residual = fitTrace.y - prediction;
    fit.rmse = sqrt(mean(fit.residual.^2));
    fit.coefficients = coefficients;
    fit.designCondition = conditionNumber;
    fit.thermalFit = design(:, 1:3) * coefficients(1:3);
    fit.sawFit = prediction - fit.thermalFit;
    fit.sawActive = saw.active;
    fit.sawFrequencyHz = saw.frequencyHz;
    fit.sawTau = NaN;
    fit.sawAtBound = false;
    if saw.active
        fit.sawTau = 10.^bestX(4);
        fit.sawAtBound = nearBound(4);
    end
    fit.fftInitialization = diagnostic;
    fit.multistartStarts = starts;
    fit.multistartTrials = trials;
end

function [prediction, coefficients, conditionNumber, design] = profilePrediction(x, trace, S, saw)
    % Evaluate the nonlinear basis and solve linear amplitudes.
    [temperature, displacement] = TwoLayerModel(trace.t, trace.Lambda, 10.^x(1:3), S);
    design = [temperature, displacement, ones(numel(trace.t), 1)];
    if saw.active
        envelope = exp(-trace.t / 10.^x(4));
        phase = 2 * pi * saw.frequencyHz * trace.t;
        design = [design, envelope .* cos(phase), envelope .* sin(phase)];
    end
    [prediction, coefficients, conditionNumber] = fitLinearComponents(design, trace.y);
end

function scale = residualScale(y)
    % Use one population standard deviation per run for conditioning
    scale = std(y(:), 1);
    if ~isfinite(scale) || scale <= 0
        error("TGS:DataScale", "The fit window has no finite signal variation.");
    end
end

function [x0, lower, upper] = physicalSearchSpace(S)
    % Validate thermal guesses/bounds and convert to log10
    p0 = S.p0(:).';
    pLower = S.pLower(:).';
    pUpper = S.pUpper(:).';
    if numel(p0) ~= 3 || numel(pLower) ~= 3 || numel(pUpper) ~= 3 || ...
            any(~isfinite([p0, pLower, pUpper])) || any([p0, pLower, pUpper] <= 0) || ...
            any(pLower >= pUpper) || any(p0 <= pLower | p0 >= pUpper)
        error("TGS:Parameters", "Positive [alpha_f,alpha_s,R] guesses must lie strictly inside bounds.");
    end
    x0 = log10(p0);
    lower = log10(pLower);
    upper = log10(pUpper);
end

function [tau0, bounds] = acousticLifetime(trace, frequencyHz)
    % Set the SAW decay-time guess and bounds [s]
    % Bounds run from half an acoustic period to twice the fitted time span
    duration = trace.t(end) - trace.t(1);
    bounds = [0.5 / frequencyHz, 2 * duration];
    if numel(bounds) ~= 2 || any(~isfinite(bounds)) || bounds(1) <= 0 || bounds(1) >= bounds(2)
        error("TGS:SAWLifetime", "Automatic SAW lifetime bounds are invalid.");
    end
    margin = 0.05 * diff(bounds);
    tau0 = min(max(duration / 5, bounds(1) + margin), bounds(2) - margin);
    if ~isscalar(tau0) || ~isfinite(tau0) || tau0 <= bounds(1) || tau0 >= bounds(2)
        error("TGS:SAWLifetime", "Automatic SAW lifetime guess must lie inside its bounds.");
    end
end

function starts = multistartPoints(primary, lower, upper, S)
    % Build reproducible starts in bounded log10 coordinates
    % Retain unique anchors up to multistartCount
    % fill remaining starts with points using a local random stream
    anchors = unique([primary; (lower + upper) / 2], "rows", "stable");
    anchors = anchors(1:min(S.multistartCount, size(anchors, 1)), :);
    count = S.multistartCount - size(anchors, 1);
    stream = RandStream("mt19937ar", "Seed", S.multistartSeed);
    unitPoints = zeros(count, numel(lower));
    for j = 1:numel(lower)
        if count > 0
            strata = randperm(stream, count).';
            unitPoints(:, j) = (strata - 1 + rand(stream, count, 1)) / count;
        end
    end
    starts = [anchors; lower + unitPoints .* (upper - lower)];
end

function [prediction, coefficients, conditionNumber] = fitLinearComponents(design, y)
    % Solve unconstrained linear amplitudes by scaled SVD
    y = y(:);
    if size(design, 1) ~= numel(y) || any(~isfinite(design), "all") || ...
            any(~isfinite(y))
        error("TGS:LinearFit", "Design and data must have paired finite rows.");
    end
    scale = vecnorm(design, 2, 1);
    if any(scale == 0)
        error("TGS:LinearFit", "A design component has zero magnitude.");
    end
    scaledDesign = design ./ scale;
    [left, singular, right] = svd(scaledDesign, "econ");
    values = diag(singular);
    keep = values > max(size(scaledDesign)) * eps(values(1));
    scaledCoefficients = right(:, keep) * ...
        ((left(:, keep)' * y) ./ values(keep));
    coefficients = scaledCoefficients ./ scale.';
    prediction = design * coefficients;
    if nargout >= 3
        if nnz(keep) < size(design, 2)
            conditionNumber = Inf;
        else
            conditionNumber = values(1) / values(end);
        end
    end
end

%% Physical Model
function [temperature, displacement] = TwoLayerModel(t, Lambda, p, S)
    % Compute film-surface thermal and displacement basis responses

    % Validate and orient model inputs.
        t = t(:);
        p = p(:).';

    % Calculate material and coupling constants.
        alpha_f = p(1);
        alpha_s = p(2);
        R = p(3);
        L = S.L;
        q = 2 * pi / Lambda;
        filmConductivity = S.Cf * alpha_f;
        substrateConductivity = S.Cs * alpha_s;
        mu_f = S.Ef / (2 * (1 + S.nuf));
        mu_s = S.Es / (2 * (1 + S.nus));
        gamma_f = S.alphathf * (1 + S.nuf) / (1 - S.nuf);
        gamma_s = S.alphaths * (1 + S.nus) / (1 - S.nus);

    % Build the symmetric angular-frequency grid.
        positiveCount = S.Nw / 2;
        omegaPositive = logspace(log10(S.wmin), log10(S.wmax), positiveCount);
        omega = [-fliplr(omegaPositive), omegaPositive];
        temperaturePositive = complex(zeros(size(omegaPositive)));
        displacementPositive = complex(zeros(size(omegaPositive)));

    % Assemble the elastic boundary-condition system.
        exp2qL = exp(2 * q * L);
        shearRatio = mu_s / mu_f;
        elasticMatrix = [ 1, 2*S.nuf-2,         1,      2-2*S.nuf,               0           0; ...
                          1, 2*S.nuf-1,        -1,      2*S.nuf-1,               0,          0; ...
                          1, q*L-3+4*S.nuf,    -exp2qL, exp2qL*(4*S.nuf-q*L-3), -1,          3-4*S.nus-q*L; ...
                          1, q*L,               exp2qL, q*L*exp2qL,             -1,         -q*L; ...
                          1, q*L-2+2*S.nuf,     exp2qL, exp2qL*(q*L+2-2*S.nuf), -shearRatio, (2-2*S.nus-q*L)*shearRatio; ...
                          1, q*L-1+2*S.nuf,    -exp2qL, exp2qL*(2*S.nuf-q*L-1), -shearRatio, (1-2*S.nus-q*L)*shearRatio];
        elasticSolver = decomposition(elasticMatrix, "lu");

    % Solve the coupled response at each positive frequency.
        for frequencyIndex = 1:positiveCount
            currentOmega = omegaPositive(frequencyIndex);
            beta_f = sqrt(q^2 + 1i * currentOmega / alpha_f);
            beta_s = sqrt(q^2 + 1i * currentOmega / alpha_s);
            r_f = 1i * currentOmega / alpha_f;
            r_s = 1i * currentOmega / alpha_s;
            conductivityRatio = substrateConductivity * beta_s / (filmConductivity * beta_f);
            interfaceFactor = 1 - conductivityRatio + substrateConductivity * beta_s * R;
            thermalDenominator = 2 * conductivityRatio + interfaceFactor * (-expm1(-2 * beta_f * L));
            commonFactor = S.Q0 / (filmConductivity * beta_f * thermalDenominator);
            interfaceDecay = exp((q - beta_f) * L);
            A_f = commonFactor * (2 * conductivityRatio + interfaceFactor);
            B_f = commonFactor * interfaceFactor * exp(-2 * beta_f * L);
            A_f_interface = commonFactor * (2 * conductivityRatio + interfaceFactor) * interfaceDecay;
            B_f_interface = commonFactor * interfaceFactor * interfaceDecay;
            A_s_interface = 2 * commonFactor * interfaceDecay;
            elasticForcing = [ gamma_f * beta_f * (A_f - B_f) / r_f; ...
                               gamma_f * q * (A_f + B_f) / r_f; ...
                               gamma_f * q * (A_f_interface + B_f_interface) / r_f - ...
                               gamma_s * q * A_s_interface / r_s; ...
                               gamma_f * beta_f * (A_f_interface - B_f_interface) / r_f - gamma_s * beta_s * A_s_interface / r_s; ...
                               gamma_f * beta_f * (A_f_interface - B_f_interface) / r_f - shearRatio * gamma_s * beta_s * A_s_interface / r_s; ...
                               gamma_f * q * (A_f_interface + B_f_interface) / r_f -  shearRatio * gamma_s * q * A_s_interface / r_s];
            elasticCoefficients = elasticSolver \ elasticForcing;
            temperaturePositive(frequencyIndex) = A_f + B_f;
            displacementPositive(frequencyIndex) = gamma_f * beta_f * (B_f - A_f) / r_f + elasticCoefficients(1) + elasticCoefficients(3);
        end

    % Reconstruct the time-domain responses
        temperatureSpectrum = [conj(fliplr(temperaturePositive)), temperaturePositive];
        displacementSpectrum = [conj(fliplr(displacementPositive)), displacementPositive];
        spectra = [temperatureSpectrum(:), displacementSpectrum(:)];
        response = real(inverseFourierPiecewiseLinear(omega, spectra, t)) / (2 * pi);
        temperature = response(:, 1);
        displacement = response(:, 2);
end

function timeSignal = inverseFourierPiecewiseLinear(omega, spectra, t)
    % Integrate spectra on a nonuniform omega grid
    omega = omega(:);
    deltaOmega = diff(omega);
    spectralSlope = diff(spectra, 1, 1) ./ deltaOmega;
    timeSignal = complex(zeros(numel(t), size(spectra, 2)));
    for timeIndex = 1:numel(t)
        currentTime = t(timeIndex);
        useSeries = abs(currentTime * deltaOmega) < 1e-4;
        integralConstant = complex(zeros(size(deltaOmega)));
        integralLinear = complex(zeros(size(deltaOmega)));
        interval = deltaOmega(~useSeries);
        phase = exp(1i * currentTime * interval);
        integralConstant(~useSeries) = expm1(1i * currentTime * interval) / (1i * currentTime);
        integralLinear(~useSeries) = phase .* (interval / (1i * currentTime) + 1 / currentTime^2) - 1 / currentTime^2;
        interval = deltaOmega(useSeries);
        integralConstant(useSeries) = interval + ...
            1i * currentTime * interval.^2 / 2 - currentTime^2 * interval.^3 / 6;
        integralLinear(useSeries) = interval.^2 / 2 + ...
            1i * currentTime * interval.^3 / 3 - currentTime^2 * interval.^4 / 8;
        timeSignal(timeIndex, :) = sum(exp(1i * currentTime * omega(1:end - 1)) .* ...
            (spectra(1:end - 1, :) .* integralConstant + spectralSlope .* integralLinear), 1);
    end
end

%% Plots and Sensitivity
function figures = plotTwoLayerTGSResults(traces, fits, S)
    figures = gobjects(numel(fits), 1);
    plotTimes = vertcat(traces.t);
    plotTimes = plotTimes(plotTimes >= S.tFit(1) & plotTimes <= S.tFit(2));
    sharedLimitsNs = [min(plotTimes), max(plotTimes)] * 1e9;
    for i = 1:numel(fits)
        fit = fits(i);
        trace = traces([traces.run] == fit.run);
        if isfield(trace, "tFull") && ~isempty(trace.tFull)
            t = trace.tFull(:);
            y = trace.yFull(:);
        else
            t = trace.t(:);
            y = trace.y(:);
        end
        figures(i) = figure("Name", sprintf("TGS run %d", trace.run));
        layout = tiledlayout(figures(i), 1, 1, "Padding", "compact");
        axisHandle = nexttile(layout);
        use = t >= S.tFit(1) & t <= S.tFit(2);
        plot(axisHandle, t(use) * 1e9, y(use), "Color", [0.65, 0.65, 0.65]);
        hold(axisHandle, "on");
        plot(axisHandle, fit.t * 1e9, fit.yFit, "r-", "LineWidth", 1.4);
        plot(axisHandle, fit.t * 1e9, fit.thermalFit, "b--", "LineWidth", 1.1);
        xline(axisHandle, fit.t(1) * 1e9, "k:");
        xlabel(axisHandle, "Time relative to pump arrival (ns)");
        ylabel(axisHandle, "Detector signal");
        legend(axisHandle, "Raw data", "Total fit", "Thermal + offset", ...
            "Fit start", "Location", "best");
         text(axisHandle, 0.98, 0.02, ...
            sprintf("RMSE = %.3g", fit.rmse), ...
            "Units", "normalized", ...
            "HorizontalAlignment", "right", ...
            "VerticalAlignment", "middle", ...
            "Interpreter", "tex", "Margin", 4);
        xlim(axisHandle, sharedLimitsNs);
        grid(axisHandle, "on");
        title(axisHandle, sprintf( ...
            "Run %d: alpha_f = %.3g, alpha_s = %.3g m^2/s; R = %.3g m^2 K/W", ...
            trace.run, fit.p(1), fit.p(2), fit.p(3)));
    end
end

function [sensitivity, figureHandle] = startingValueSensitivity(traces, S, baselineResults, initializationTraces)
    % Refit after perturbing one thermal start at a time
    % Change p0 by low/high factors and rerun 
    % Store per-run percent changes and a median summary in table UserData
    % Tests dependence on starting values

    parameters = ["AlphaF", "AlphaS", "R"];
    factors = S.sensitivityStartFactors(:).';
    if numel(factors) ~= 2 || any(~isfinite(factors)) || ...
            factors(1) <= 0 || factors(1) >= 1 || factors(2) <= 1
        error("TGS:Sensitivity", ...
            "sensitivityStartFactors must satisfy 0 < low < 1 < high.");
    end
    starts = S.p0(:) .* factors;
    if any(starts <= S.pLower(:) | starts >= S.pUpper(:), "all")
        error("TGS:Sensitivity", ...
            "Perturbed starts must lie strictly inside the physical bounds.");
    end

    baselineValues = baselineResults{:, cellstr(parameters)};
    rows = cell(18, 1);
    row = 0;
    for changed = 1:3
        for factorIndex = 1:2
            scenario = S;
            scenario.p0(changed) = S.p0(changed) * factors(factorIndex);
            [scenarioResults, ~] = fitTwoLayerTGSRuns( ...
                traces, scenario, initializationTraces);
            [found, order] = ismember(baselineResults.Run, scenarioResults.Run);
            if ~all(found)
                error("TGS:Sensitivity", ...
                    "Sensitivity and baseline run identifiers differ.");
            end
            values = scenarioResults{order, cellstr(parameters)};
            for response = 1:3
                row = row + 1;
                count = height(baselineResults);
                changePercent = 100 * (values(:, response) ./ ...
                    baselineValues(:, response) - 1);
                rows{row} = table(baselineResults.Run, ...
                    repmat(parameters(changed), count, 1), ...
                    repmat(parameters(response), count, 1), ...
                    repmat(factors(factorIndex), count, 1), ...
                    baselineValues(:, response), values(:, response), ...
                    changePercent, scenarioResults.RMSE(order), ...
                    scenarioResults.Converged(order), ...
                    'VariableNames', {'Run', 'PerturbedParameter', ...
                    'ResponseParameter', 'StartFactor', 'BaselineFit', ...
                    'PerturbedFit', 'FitChangePercent', 'RMSE', 'Converged'});
            end
        end
    end
    sensitivity = vertcat(rows{:});

    summaryRows = cell(18, 1);
    row = 0;
    for response = 1:3
        for changed = 1:3
            for factorIndex = 1:2
                use = sensitivity.ResponseParameter == parameters(response) & ...
                    sensitivity.PerturbedParameter == parameters(changed) & ...
                    sensitivity.StartFactor == factors(factorIndex);
                row = row + 1;
                summaryRows{row} = table(parameters(changed), ...
                    parameters(response), factors(factorIndex), ...
                    median(sensitivity.FitChangePercent(use)), ...
                    'VariableNames', {'PerturbedParameter', ...
                    'ResponseParameter', 'StartFactor', ...
                    'MedianFitChangePercent'});
            end
        end
    end
    summary = vertcat(summaryRows{:});
    sensitivity.Properties.UserData.Summary = summary;

    figureHandle = gobjects(0);
    if ~S.showPlots
        return
    end

    figureHandle = figure("Name", "Starting-value sensitivity tornado");
    layout = tiledlayout(figureHandle, 1, 3, "TileSpacing", "compact", ...
        "Padding", "compact");
    responseLabels = ["Film diffusivity", "Substrate diffusivity", ...
        "Interface resistance"];
    parameterLabels = ["alpha_f", "alpha_s", "R"];
    legendBars = gobjects(0);
    legendAxis = gobjects(0);
    for response = 1:3
        effects = zeros(3, 2);
        for changed = 1:3
            for factorIndex = 1:2
                use = summary.ResponseParameter == parameters(response) & ...
                    summary.PerturbedParameter == parameters(changed) & ...
                    summary.StartFactor == factors(factorIndex);
                effects(changed, factorIndex) = ...
                    summary.MedianFitChangePercent(use);
            end
        end
        [~, order] = sort(abs(diff(effects, 1, 2)), "ascend");
        axisHandle = nexttile(layout);
        bars = barh(axisHandle, effects(order, :), "grouped");
        bars(1).FaceColor = [0.25, 0.48, 0.75];
        bars(2).FaceColor = [0.85, 0.40, 0.25];
        xline(axisHandle, 0, "k-");
        axisHandle.YTick = 1:3;
        axisHandle.YTickLabel = parameterLabels(order);
        limit = max(abs(effects), [], "all");
        if limit == 0
            limit = 1;
        end
        axisHandle.XLim = 1.1 * [-limit, limit];
        xlabel(axisHandle, "Median fitted change (%)");
        title(axisHandle, responseLabels(response));
        grid(axisHandle, "on");
        if response == 1
            legendBars = bars;
            legendAxis = axisHandle;
        end
    end
    title(layout, "Starting-value sensitivity");
    legendHandle = legend(legendAxis, legendBars, ...
        compose("Start x %.2g", factors), "Orientation", "horizontal");
    legendHandle.Layout.Tile = "south";
end
