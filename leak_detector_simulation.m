%% Моделирование концентрационных полей при утечке в замкнутом помещении
% Расширенная версия: сценарный анализ расхода утечки, вентиляции,
% воздушных потоков и случайных возмущений. Для каждой ситуации строятся
% отдельные графики, чтобы наглядно показать отличие от базового режима.

clear; clc; close all;

%% 1) Геометрия помещения и расчетная сетка
params.Lx = 10;                  % длина помещения, м
params.Ly = 6;                   % ширина помещения, м
params.Nx = 140;                 % число узлов по оси x
params.Ny = 84;                  % число узлов по оси y
params.dx = params.Lx/(params.Nx-1);
params.dy = params.Ly/(params.Ny-1);
[params.x, params.y] = meshgrid(linspace(0, params.Lx, params.Nx), ...
                                linspace(0, params.Ly, params.Ny));

%% 2) Общие физические параметры и критерии безопасности
params.D = 1.4e-2;               % эффективный коэффициент диффузии, м^2/с
params.Tend = 240;               % длительность моделирования, с
params.xLeak = 2.4;
params.yLeak = 1.7;
params.sigma = 0.22;
params.dangerThreshold = 0.050;  % порог опасной концентрации, усл. ед.
params.safeThreshold = 0.018;    % порог безопасного режима, усл. ед.
params.safeHoldTime = 20;        % удержание ниже safeThreshold, с
params.maxExpectedSpeed = 0.12;  % запас для устойчивого шага при порывах, м/с

%% 3) Источник утечки и траектория течеискателя
src = exp(-((params.x-params.xLeak).^2 + (params.y-params.yLeak).^2)/(2*params.sigma^2));
params.src = src / sum(src(:));
[params.pathX, params.pathY] = buildScanPath(params.Lx, params.Ly, params.Tend);

%% 4) Базовый сценарий
baseScenario = makeScenario('Базовый режим', 8.0, 0.02, 0.00, 0.000, inf, 'uniform', struct([]));
baseResult = runScenario(baseScenario, params);

%% 5) Ситуация 1: рост расхода утечки Q
leakScenarios = [ ...
    makeScenario('Малый расход Q = 4', 4.0, 0.02, 0.00, 0.000, inf, 'uniform', struct([])), ...
    baseScenario, ...
    makeScenario('Большой расход Q = 12', 12.0, 0.02, 0.00, 0.000, inf, 'uniform', struct([]))];
leakResults = runScenarioSet(leakScenarios, params);
plotLeakRateComparison(leakResults, params);

%% 6) Ситуация 2: повышение производительности вентиляции
% Для анализа выхода на безопасный режим задаем прекращение утечки после 120 с,
% а затем сравниваем, как быстро разные уровни вентиляции удаляют вещество.
ventScenarios = [ ...
    makeScenario('Слабая вентиляция k = 0.003 1/с', 8.0, 0.02, 0.00, 0.003, 120, 'uniform', struct([])), ...
    makeScenario('Средняя вентиляция k = 0.010 1/с', 8.0, 0.02, 0.00, 0.010, 120, 'uniform', struct([])), ...
    makeScenario('Сильная вентиляция k = 0.025 1/с', 8.0, 0.02, 0.00, 0.025, 120, 'uniform', struct([]))];
ventResults = runScenarioSet(ventScenarios, params);
plotVentilationComparison(ventResults, params);

%% 7) Ситуация 3: асимметрия от воздушных потоков
flowScenarios = [ ...
    makeScenario('Без направленного потока', 8.0, 0.00, 0.00, 0.000, inf, 'uniform', struct([])), ...
    makeScenario('Равномерный поток вправо', 8.0, 0.04, 0.00, 0.000, inf, 'uniform', struct([])), ...
    makeScenario('Диагональный поток', 8.0, 0.03, 0.025, 0.000, inf, 'uniform', struct([])), ...
    makeScenario('Локальный сквозняк', 8.0, 0.02, 0.00, 0.000, inf, 'localDraft', struct([]))];
flowResults = runScenarioSet(flowScenarios, params);
plotAirflowComparison(flowResults, params);

%% 8) Ситуация 4: случайные возмущения
rng(7); % воспроизводимость набора возмущений
noDisturbance = makeScenario('Без возмущений', 8.0, 0.02, 0.00, 0.006, inf, 'uniform', struct([]));
disturbanceScenarios = [ ...
    noDisturbance, ...
    makeScenario('Открывание двери', 8.0, 0.02, 0.00, 0.006, inf, 'uniform', makeDoorDisturbance()), ...
    makeScenario('Изменение направления потока', 8.0, 0.02, 0.00, 0.006, inf, 'uniform', makeFlowDirectionDisturbance()), ...
    makeScenario('Локальный сквозняк', 8.0, 0.02, 0.00, 0.006, inf, 'uniform', makeLocalDraftDisturbance()), ...
    makeScenario('Изменение температуры', 8.0, 0.02, 0.00, 0.006, inf, 'uniform', makeTemperatureDisturbance()), ...
    makeScenario('Изменение давления', 8.0, 0.02, 0.00, 0.006, inf, 'uniform', makePressureDisturbance())];
disturbanceResults = runScenarioSet(disturbanceScenarios, params);
plotDisturbanceComparison(disturbanceResults, params);

%% 9) Базовая визуализация течеискателя
plotDetectorResponse(baseResult, params);

fprintf('\nСводка сценариев:\n');
printSummary([leakResults, ventResults, flowResults, disturbanceResults], params);

%% Локальные функции
function scenario = makeScenario(name, Q, uX, uY, ventRate, leakStopTime, flowMode, disturbances)
    scenario.name = name;
    scenario.Q = Q;
    scenario.uX = uX;
    scenario.uY = uY;
    scenario.ventRate = ventRate;
    scenario.leakStopTime = leakStopTime;
    scenario.flowMode = flowMode;
    scenario.disturbances = disturbances;
end

function [pathX, pathY] = buildScanPath(Lx, Ly, Tend)
    scanLines = 8;
    yLines = linspace(0.5, Ly-0.5, scanLines);
    pathX = [];
    pathY = [];
    for k = 1:scanLines
        if mod(k,2)==1
            pathX = [pathX, linspace(0.4, Lx-0.4, 220)]; %#ok<AGROW>
        else
            pathX = [pathX, linspace(Lx-0.4, 0.4, 220)]; %#ok<AGROW>
        end
        pathY = [pathY, yLines(k)*ones(1,220)]; %#ok<AGROW>
    end
    pathT = linspace(0, Tend, numel(pathX));
    sampleT = linspace(0, Tend, ceil(Tend/0.1));
    pathX = interp1(pathT, pathX, sampleT, 'linear', 'extrap');
    pathY = interp1(pathT, pathY, sampleT, 'linear', 'extrap');
end

function results = runScenarioSet(scenarios, params)
    results = repmat(runScenario(scenarios(1), params), 1, numel(scenarios));
    for s = 1:numel(scenarios)
        results(s) = runScenario(scenarios(s), params);
    end
end

function result = runScenario(scenario, params)
    dtDiff = 0.25*min(params.dx^2, params.dy^2)/params.D;
    dtAdv = 0.45/max(params.maxExpectedSpeed/min(params.dx, params.dy), 1e-12);
    dt = min(dtDiff, dtAdv);
    Nt = ceil(params.Tend/dt);
    time = (0:Nt-1)*dt;

    C = zeros(params.Ny, params.Nx);
    meanC = zeros(1, Nt);
    maxC = zeros(1, Nt);
    dangerArea = zeros(1, Nt);
    sensorSignal = zeros(1, Nt);
    sensorSampleT = linspace(0, params.Tend, numel(params.pathX));
    sensorX = interp1(sensorSampleT, params.pathX, time, 'linear', 'extrap');
    sensorY = interp1(sensorSampleT, params.pathY, time, 'linear', 'extrap');

    for n = 1:Nt
        t = time(n);
        current = applyDisturbances(scenario, params, t);
        activeQ = current.Q;
        if t > current.leakStopTime
            activeQ = 0;
        end

        dCdx = zeros(size(C)); dCdy = zeros(size(C));
        d2Cdx2 = zeros(size(C)); d2Cdy2 = zeros(size(C));

        dCdx(:,2:end-1) = (C(:,3:end)-C(:,1:end-2))/(2*params.dx);
        dCdy(2:end-1,:) = (C(3:end,:)-C(1:end-2,:))/(2*params.dy);
        d2Cdx2(:,2:end-1) = (C(:,3:end)-2*C(:,2:end-1)+C(:,1:end-2))/params.dx^2;
        d2Cdy2(2:end-1,:) = (C(3:end,:)-2*C(2:end-1,:)+C(1:end-2,:))/params.dy^2;

        C = C + dt*(current.D*(d2Cdx2+d2Cdy2) ...
            - current.Ux.*dCdx - current.Uy.*dCdy ...
            + activeQ*params.src - current.ventRate*C);
        C = max(C, 0);

        C(:,1) = C(:,2);
        C(:,end) = C(:,end-1);
        C(1,:) = C(2,:);
        C(end,:) = C(end-1,:);

        sensorSignal(n) = interpolateSensor(C, sensorX(n), sensorY(n), params);
        meanC(n) = mean(C(:));
        maxC(n) = max(C(:));
        dangerArea(n) = sum(C(:) > params.dangerThreshold)*params.dx*params.dy;
    end

    finalFlow = applyDisturbances(scenario, params, params.Tend);
    result.name = scenario.name;
    result.scenario = scenario;
    result.time = time;
    result.dt = dt;
    result.C = C;
    result.Ux = finalFlow.Ux;
    result.Uy = finalFlow.Uy;
    result.meanC = meanC;
    result.maxC = maxC;
    result.dangerArea = dangerArea;
    result.sensorSignal = sensorSignal;
    result.sensorX = sensorX;
    result.sensorY = sensorY;
    result.firstDangerTime = firstTimeAbove(time, dangerArea, 0);
    result.safeModeTime = firstSustainedSafeTime(time, maxC, params.safeThreshold, params.safeHoldTime, scenario.leakStopTime);
end

function current = applyDisturbances(scenario, params, t)
    current = scenario;
    current.D = params.D;
    [current.Ux, current.Uy] = buildVelocityField(scenario.flowMode, scenario.uX, scenario.uY, params);

    for k = 1:numel(scenario.disturbances)
        d = scenario.disturbances(k);
        if t < d.startTime || t > d.startTime + d.duration
            continue;
        end
        switch d.type
            case 'doorOpen'
                current.ventRate = current.ventRate + d.extraVentRate;
                doorMask = params.x > params.Lx-1.1;
                current.Ux(doorMask) = current.Ux(doorMask) + d.extraUx;
            case 'flowDirectionChange'
                current.Ux = d.multiplierX*current.Ux;
                current.Uy = current.Uy + d.extraUy;
            case 'localDraft'
                mask = params.x > d.xMin & params.x < d.xMax & params.y > d.yMin & params.y < d.yMax;
                current.Ux(mask) = current.Ux(mask) + d.extraUx;
                current.Uy(mask) = current.Uy(mask) + d.extraUy;
            case 'temperatureChange'
                current.D = current.D*d.diffusionMultiplier;
            case 'pressureChange'
                current.Q = current.Q*d.sourceMultiplier;
                current.ventRate = max(0, current.ventRate + d.extraVentRate);
        end
    end
end

function [Ux, Uy] = buildVelocityField(flowMode, uX, uY, params)
    Ux = uX*ones(params.Ny, params.Nx);
    Uy = uY*ones(params.Ny, params.Nx);
    switch flowMode
        case 'uniform'
            return;
        case 'localDraft'
            draftMask = params.x > 5.6 & params.x < 8.8 & params.y > 2.0 & params.y < 3.4;
            Ux(draftMask) = Ux(draftMask) + 0.065;
            Uy(draftMask) = Uy(draftMask) - 0.020;
        otherwise
            error('Неизвестный режим потока: %s', flowMode);
    end
end

function value = interpolateSensor(C, sx, sy, params)
    i = min(max(floor(sx/params.dx)+1,1),params.Nx-1);
    j = min(max(floor(sy/params.dy)+1,1),params.Ny-1);
    tx = (sx - (i-1)*params.dx)/params.dx;
    ty = (sy - (j-1)*params.dy)/params.dy;
    c00 = C(j,i);   c10 = C(j,i+1);
    c01 = C(j+1,i); c11 = C(j+1,i+1);
    value = (1-tx)*(1-ty)*c00 + tx*(1-ty)*c10 + (1-tx)*ty*c01 + tx*ty*c11;
end

function t = firstTimeAbove(time, values, threshold)
    idx = find(values > threshold, 1, 'first');
    if isempty(idx)
        t = NaN;
    else
        t = time(idx);
    end
end

function t = firstSustainedSafeTime(time, values, threshold, holdTime, startTime)
    t = NaN;
    if numel(time) < 2 || ~isfinite(startTime)
        return;
    end
    dt = time(2)-time(1);
    holdSteps = max(1, ceil(holdTime/dt));
    below = values < threshold;
    startIdx = find(time >= startTime, 1, 'first');
    if isempty(startIdx)
        return;
    end
    for k = startIdx:numel(values)-holdSteps+1
        if all(below(k:k+holdSteps-1))
            t = time(k);
            return;
        end
    end
end

function d = makeDoorDisturbance()
    d = struct('type','doorOpen','startTime',70,'duration',45, ...
        'extraVentRate',0.018,'extraUx',0.055);
end

function d = makeFlowDirectionDisturbance()
    d = struct('type','flowDirectionChange','startTime',85,'duration',65, ...
        'multiplierX',-1.0,'extraUy',0.040);
end

function d = makeLocalDraftDisturbance()
    d = struct('type','localDraft','startTime',95,'duration',50, ...
        'xMin',5.2,'xMax',8.9,'yMin',2.0,'yMax',3.6,'extraUx',0.075,'extraUy',-0.030);
end

function d = makeTemperatureDisturbance()
    d = struct('type','temperatureChange','startTime',100,'duration',70, ...
        'diffusionMultiplier',1.75);
end

function d = makePressureDisturbance()
    d = struct('type','pressureChange','startTime',110,'duration',55, ...
        'sourceMultiplier',1.55,'extraVentRate',-0.002);
end

function plotDetectorResponse(result, params)
    figure('Name','Базовый отклик течеискателя','Color','w','Position',[80 80 1200 520]);
    subplot(1,2,1);
    imagesc([0 params.Lx],[0 params.Ly],result.C);
    set(gca,'YDir','normal'); axis equal tight; hold on;
    plot(params.pathX, params.pathY, 'w--', 'LineWidth', 1.0);
    plot(result.sensorX(end), result.sensorY(end), 'wo', 'MarkerFaceColor','k', 'MarkerSize', 6);
    plot(params.xLeak, params.yLeak, 'rp', 'MarkerFaceColor','r', 'MarkerSize', 14);
    colorbar; xlabel('x, м'); ylabel('y, м');
    title('Итоговое поле концентрации');
    legend({'Траектория прибора','Текущее положение','Точка утечки'}, 'Location','northoutside');

    subplot(1,2,2);
    plot(result.time, result.sensorSignal, 'b-', 'LineWidth', 1.5); grid on;
    xlabel('Время, с'); ylabel('Концентрация, усл. ед.');
    title('Отклик течеискателя во время сканирования');
    sgtitle('Базовый режим локализации утечки');
end

function plotLeakRateComparison(results, params)
    figure('Name','Ситуация 1: расход утечки','Color','w','Position',[70 70 1250 760]);
    subplot(2,2,1); hold on; grid on;
    for k = 1:numel(results)
        plot(results(k).time, results(k).dangerArea, 'LineWidth', 1.5);
    end
    xlabel('Время, с'); ylabel('Площадь опасной зоны, м^2');
    title('Чем больше Q, тем быстрее растет опасная зона');
    legend({results.name}, 'Location','northwest');

    subplot(2,2,2); hold on; grid on;
    for k = 1:numel(results)
        plot(results(k).time, results(k).maxC, 'LineWidth', 1.5);
    end
    yline(params.dangerThreshold, 'r--', 'Порог опасности');
    xlabel('Время, с'); ylabel('Максимальная концентрация');
    title('Максимальная концентрация при разных расходах');

    subplot(2,2,3);
    times = [results.firstDangerTime];
    bar(times); grid on;
    set(gca,'XTick',1:numel(results),'XTickLabel',{results.name},'XTickLabelRotation',20);
    ylabel('Время, с'); title('Начало формирования опасной зоны');

    subplot(2,2,4); hold on; axis equal tight;
    imagesc([0 params.Lx],[0 params.Ly],results(end).C-results(1).C);
    set(gca,'YDir','normal'); colorbar;
    plot(params.xLeak, params.yLeak, 'rp', 'MarkerFaceColor','r', 'MarkerSize', 12);
    xlabel('x, м'); ylabel('y, м');
    title('Разница полей: большой Q минус малый Q');
    sgtitle('Отдельный график влияния расхода утечки');
end

function plotVentilationComparison(results, params)
    figure('Name','Ситуация 2: вентиляция','Color','w','Position',[90 90 1250 760]);
    subplot(2,2,1); hold on; grid on;
    for k = 1:numel(results)
        plot(results(k).time, results(k).meanC, 'LineWidth', 1.5);
    end
    xline(120, 'k--', 'Остановка утечки');
    xlabel('Время, с'); ylabel('Средняя концентрация');
    title('Вентиляция снижает среднюю концентрацию');
    legend({results.name}, 'Location','northwest');

    subplot(2,2,2); hold on; grid on;
    for k = 1:numel(results)
        plot(results(k).time, results(k).maxC, 'LineWidth', 1.5);
    end
    yline(params.safeThreshold, 'g--', 'Безопасный порог');
    xline(120, 'k--', 'Остановка утечки');
    xlabel('Время, с'); ylabel('Максимальная концентрация');
    title('Время выхода на безопасный режим');

    subplot(2,2,3);
    safeTimes = [results.safeModeTime];
    bar(safeTimes); grid on;
    set(gca,'XTick',1:numel(results),'XTickLabel',{results.name},'XTickLabelRotation',20);
    ylabel('Время, с'); title('Чем выше вентиляция, тем меньше время очистки');

    subplot(2,2,4); hold on; grid on;
    reference = results(1).meanC;
    for k = 2:numel(results)
        plot(results(k).time, results(k).meanC-reference, 'LineWidth', 1.5);
    end
    xlabel('Время, с'); ylabel('\Delta средней концентрации');
    title('Разница относительно слабой вентиляции');
    legend({results(2:end).name}, 'Location','southwest');
    sgtitle('Отдельный график влияния производительности вентиляции');
end

function plotAirflowComparison(results, params)
    figure('Name','Ситуация 3: воздушные потоки','Color','w','Position',[110 110 1300 780]);
    baseC = results(1).C;
    for k = 1:numel(results)
        subplot(2,numel(results),k);
        imagesc([0 params.Lx],[0 params.Ly],results(k).C);
        set(gca,'YDir','normal'); axis equal tight; hold on;
        quiver(params.x(1:8:end,1:8:end), params.y(1:8:end,1:8:end), ...
               results(k).Ux(1:8:end,1:8:end), results(k).Uy(1:8:end,1:8:end), 'k');
        plot(params.xLeak, params.yLeak, 'rp', 'MarkerFaceColor','r', 'MarkerSize', 10);
        colorbar; title(results(k).name); xlabel('x, м'); ylabel('y, м');

        subplot(2,numel(results),numel(results)+k);
        imagesc([0 params.Lx],[0 params.Ly],results(k).C-baseC);
        set(gca,'YDir','normal'); axis equal tight; colorbar;
        xlabel('x, м'); ylabel('y, м');
        title('\DeltaC относительно режима без потока');
    end
    sgtitle('Отдельный график асимметрии концентрационного поля от воздушных потоков');
end

function plotDisturbanceComparison(results, params)
    baseline = results(1);
    figure('Name','Ситуация 4: случайные возмущения','Color','w','Position',[130 130 1320 860]);
    disturbanceCount = numel(results)-1;
    for k = 2:numel(results)
        row = k-1;
        subplot(disturbanceCount,3,3*(row-1)+1); hold on; grid on;
        plot(baseline.time, baseline.meanC, 'k--', 'LineWidth', 1.0);
        plot(results(k).time, results(k).meanC, 'LineWidth', 1.3);
        xlabel('t, с'); ylabel('Средняя C'); title(results(k).name);
        if row == 1
            legend({'Без возмущений','С возмущением'}, 'Location','northwest');
        end

        subplot(disturbanceCount,3,3*(row-1)+2); hold on; grid on;
        plot(results(k).time, results(k).maxC-baseline.maxC, 'LineWidth', 1.3);
        xlabel('t, с'); ylabel('\Delta max(C)'); title('Разница максимума');

        subplot(disturbanceCount,3,3*(row-1)+3);
        imagesc([0 params.Lx],[0 params.Ly],results(k).C-baseline.C);
        set(gca,'YDir','normal'); axis equal tight; colorbar;
        xlabel('x, м'); ylabel('y, м'); title('Разница итогового поля');
    end
    sgtitle('Отдельные графики влияния случайных возмущений');
end

function printSummary(results, params)
    seen = containers.Map();
    for k = 1:numel(results)
        key = results(k).name;
        if isKey(seen, key)
            continue;
        end
        seen(key) = true;
        fprintf('  %-38s | maxC=%7.4f | meanC=%7.4f | опасная зона с t=%7.2f c | безопасный режим t=%7.2f c\n', ...
            results(k).name, max(results(k).maxC), results(k).meanC(end), ...
            results(k).firstDangerTime, results(k).safeModeTime);
    end
    fprintf('  Пороги: dangerThreshold=%.3f, safeThreshold=%.3f, safeHoldTime=%.1f c\n', ...
        params.dangerThreshold, params.safeThreshold, params.safeHoldTime);
end
