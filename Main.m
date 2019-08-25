clear
close all

% Тип шаблона спокойного Солнца
% 'M' - модельный
% 'T' - запись

templateType = 'T';

%% read .fits file and header
fileName = '20150715_121933_sun0_out_edit.fits';
% data
data = fitsread(fileName);
% Нужно обрезать моменты включения генератора шума
cutDataLeft = 200;
cutDataRight = 200;
data = data(:,1+cutDataLeft:end-cutDataRight,:);
% deliting bad frequencies
badFreq = 7;
%data = data(:,:,1+badFreq:end);
% header
info = fitsinfo(fileName);

% чтение значений частот из заголовка
endSuccess = 0;
commentSuccess = 0;
freqNum = 1;
counter = 1;
try
    while ~endSuccess
        if isequal(info.PrimaryData.Keywords{counter, 1}, 'COMMENT')
            commentSuccess = 1;
        end

        if isequal(info.PrimaryData.Keywords{counter+1, 1}, 'END')
            endSuccess = 1;
        end

        counter = counter + 1;

        if commentSuccess == 1 && endSuccess == 0
            freq(freqNum) = info.PrimaryData.Keywords{counter, 2};

            freqNum = freqNum + 1;
        end
    end

    % deleting bad frequencies
    freq = freq(1+badFreq:end);
catch exeption
    freq = 1:size(data,3)-badFreq;
end

% Чтение констант из заголовка fits. файла
CDELT1Success = 0;
CRPIX1Success = 0;
SOLAR_RSuccess = 0;
SOL_RASuccess = 0;
SOL_DECSuccess = 0;
i = 1;

while (i <= size(info.PrimaryData.Keywords, 1)) && ~(CDELT1Success && CRPIX1Success && SOLAR_RSuccess && SOL_RASuccess && SOL_DECSuccess)
    if isequal(info.PrimaryData.Keywords{i,1}, 'CDELT1')
        CDELT1 = info.PrimaryData.Keywords{i,2};
        CDELT1Success = 1;
    end
    
    if isequal(info.PrimaryData.Keywords{i,1}, 'CRPIX1')
        CRPIX1 = info.PrimaryData.Keywords{i,2};
        CRPIX1Success = 1;
    end
    
    if isequal(info.PrimaryData.Keywords{i,1}, 'SOLAR_R')
        SOLAR_R = info.PrimaryData.Keywords{i,2};
        SOLAR_RSuccess = 1;
    end
    
    if isequal(info.PrimaryData.Keywords{i,1}, 'SOL_RA')
        SOL_RA = info.PrimaryData.Keywords{i,2};
        SOL_RASuccess = 1;
    end
    
    if isequal(info.PrimaryData.Keywords{i,1}, 'SOL_DEC')
        SOL_DEC = info.PrimaryData.Keywords{i,2};
        SOL_DECSuccess = 1;
    end
    
    i = i + 1;

end

% Солнечный радиус из заголовка fits. файла (в пикселях)
R = 1.1*SOLAR_R / CDELT1;

if templateType == 'M'
    % Массивы значений частоты f и горизонтальной полуширины диаграммы 
    % направленности D
    % Массив взят из программы WorkScan
    f = [  0.985;   1.015;   1.045;   1.670;   1.760;   1.860;   1.950; ...
           2.050;   2.150;   2.270;   2.610;   2.720;   2.830;   2.950; ...
           3.080;   3.210;   3.350;   3.480;   3.670;   3.950;   4.270;   4.600; ...
           4.950;   5.700;   6.080;   6.500;   6.950;   7.350;   7.830; ...
           8.400;   8.750;   9.350;   9.800;  10.350;  10.950;  11.250; ...
          12.950;  13.400;  14.250;  14.750;  15.650;  16.400         ];
    D = [237.410; 235.200; 233.000; 188.000; 181.810; 174.940; 168.760; ...
         162.240; 156.070; 148.660; 128.680; 122.880; 117.090; 110.770; ...
         104.660;  99.010;  92.930;  87.280;  80.500;  70.770;  61.770;  53.580; ...
          46.550;  37.200;  33.250;  31.330;  29.270;  28.560;  27.910; ...
          27.420;  27.190;  26.730;  26.360;  25.740;  24.970;  24.500; ...
          21.580;  20.790;  19.330;  18.510;  17.290;  16.360         ] / CDELT1;

    % Эта часть находит массив горизонтальной полуширины диаграммы
    % направленности Dfreq соответсвующий массиву частот считанному из
    % fits. файла
    Dfreq(1:length(freq)) = 0;
    j = 1;
    for i = 1:length(f)
        if j <= length(freq) && freq(j) == f(i)
            Dfreq(j) = D(i);

            j = j + 1;
        end
    end

    % Получение массива вертикальной полуширины диаграммы направленности 
    % DVertFreq
    % Интерполяция взята сайта spbf.sao
    DVertFreq = 225./freq.*60;

    % Эта часть строит модельное изображение Солнца в виде квадратной
    % матрицы, которая имеет размер равный двум солнечным радиусам
    % Солнце имеет форму круга с одинаковой интенсивностью равной 1
    sun(1:4*R,1:4*R) = 0;

    % Центр круга совпадает с центром матрицы
    for i = 1:4*R
        for j = 1:4*R
            if sqrt( (i-2*R+1)^2 + (j-2*R+1)^2 ) <= R
                sun(i,j) = 1;
            end
        end
    end

    % 
    t = -2*R+1:2*R;
    % Инициализация cвертки Солнца с вертикальной диаграммой направленности
    sunShape(1:4*R,1:length(freq)) = 0;

    % Свертка Солнца с вертикальной диаграммой направленности
    for j = 1:length(freq)
        vertAntennaPattern = @(x) exp( - x.^2 / ( 2 * (DVertFreq(j)/2.355)^2 ) );

        for i = 1:4*R
            sunShape(i,j) = sum( sun(:,i)' .* vertAntennaPattern(t) );
        end
    end

    % Свертка Солнца с горизонтальной диаграммой направленности
    x = 1:4096;
    convolution(1:length(freq),1:length(x)) = 0;

    for freqNum = 1:length(freq)
        antennaPattern = @(x) exp( - (x - length(x)/2).^2 / ( 2 * (Dfreq(freqNum)/2.355)^2 ) );

        for j = 1:4*R
            convolution(freqNum,fix(j+length(x)/2 - 2*R)) = trapz(sunShape(:,freqNum)' .* antennaPattern(x(j)-t));
        end
    end
elseif templateType == 'T'
    convolution = fitsread('QuietSun201902\template.fits')';
end

%%
% Число итераций, число вписываемых гауссиан
numGauss = 8;
% Массив новых центров Солнца, на каждой частоте погоняется независимо
CRPIXfreq(1:length(freq)) = 0;
% Параметр гребневой регуляризации
% Возможно, использоваться не будет (незначительный вклад)
alpha = 0;
% Инициализация массивов с параметрами гауссиан
% Положение максимума
maxGauss(1:length(freq), 1:numGauss) = 0;
% Значение максимума
AGauss(1:length(freq), 1:numGauss) = 0;
% Полуширина
DGauss(1:length(freq), 1:numGauss) = 0;

% Создание видеопотока
v = VideoWriter('peaks.avi');
open(v)

for freqNum = 1:length(freq)
    [CRPIXfreq(freqNum), maxGauss(freqNum,:), AGauss(freqNum,:), DGauss(freqNum,:), Discr] = ...
        SunCentering(data, CRPIX1-cutDataLeft, freqNum+badFreq, numGauss, alpha, convolution(freqNum,:), 'None', v);
end

close(v)

%%
% Сортировка гауссиан
% Спектр без обработки
figure
hold on
for i = 1:numGauss
    plot(freq, maxGauss(:,i))
end
figure
hold on
for i = 1:numGauss
    plot(freq, AGauss(:,i))
    ylim([0 20000])
end
figure
hold on
for i = 1:numGauss
    plot(freq, DGauss(:,i))
end
% выбор нулей
j = 1;
for i = 1:length(freq)
    if maxGauss(i, 1) ~= 0
        freq1(j) = freq(i);
        maxGauss1(j,1:numGauss) = maxGauss(i,:);
        AGauss1(j,1:numGauss) = AGauss(i,:);
        DGauss1(j,1:numGauss) = DGauss(i,:);
        j = j + 1;
    end
end
figure
hold on
for i = 1:numGauss
    plot(freq1, maxGauss1(:,i))
end
figure
hold on
for i = 1:numGauss
    plot(freq1, AGauss1(:,i))
    ylim([0 20000])
end
figure
hold on
for i = 1:numGauss
    plot(freq1, DGauss1(:,i))
end
% сортировка по положению максимума
freq2 = freq1;
AGauss2(1:length(freq1), 1:numGauss) = 0;
DGauss2(1:length(freq1), 1:numGauss) = 0;
maxGauss2(1:length(freq1), 1:numGauss) = 0;
AGauss2(1, 1:numGauss) = AGauss1(1, :);
DGauss2(1, 1:numGauss) = DGauss1(1, :);
maxGauss2(1, 1:numGauss) = maxGauss1(1, :);
% максимально допустимое значение отклонения максимума междусоседними
% частотами, если отклонение меньше допустимого не найдено, то значение
% записывается как новый источник
trig = 10;
% размер нового массива
size2 = numGauss;

for i = 2:length(freq1)
    size1 = size2;
    
    for j = 1:numGauss
        min = trig;
        minInd = 0;
        
        for l = 1:i-1
            for k = 1:size1
                if abs(maxGauss2(l, k) - maxGauss1(i, j)) < min
                    min = abs(maxGauss2(l, k) - maxGauss1(i, j));
                    minInd = k;
                end
            end
        end
        
        if min == trig
            size2 = size2 + 1;
            AGauss2(i, size2) = AGauss1(i, j);
            DGauss2(i, size2) = DGauss1(i, j);
            maxGauss2(i, size2) = maxGauss1(i, j);
        else
            AGauss2(i, minInd) = AGauss1(i, j);
            DGauss2(i, minInd) = DGauss1(i, j);
            maxGauss2(i, minInd) = maxGauss1(i, j);
        end
    end
    
    %for j = 1:size2
    %    if AGauss2(i,j) == 0
    %        AGauss2(i,j) = AGauss2(i-1,j);
    %        DGauss2(i,j) = DGauss2(i-1,j);
    %        maxGauss2(i,j) = maxGauss2(i-1,j);
    %    end
    %end
end

for i = 1:size(AGauss2,2)
    k = 1;
    
    for j = 1:size(AGauss2,1)
        if AGauss2(j, i) ~= 0
            AGauss3{i}(k) = AGauss2(j, i);
            DGauss3{i}(k) = DGauss2(j, i);
            maxGauss3{i}(k) = maxGauss2(j, i);
            freq3{i}(k) = j;
            
            k = k + 1;
        end
    end
end

figure
hold on

% Ограничение на длину спектра (колличество точек)
lengthTrig = 15;

for i = 1:size1
    if length(freq3{i}) > lengthTrig
        plot(freq3{i}, maxGauss3{i})
    end
end
xlabel('Номер частоты')
ylabel('Положение максимума, pix')

figure
hold on
for i = 1:size1
    if length(freq3{i}) > lengthTrig
        plot(freq3{i}, AGauss3{i})
    end
end
xlabel('Номер частоты')
ylabel('Значение максимума, К')

figure
hold on
for i = 1:size1
    if length(freq3{i}) > lengthTrig
        plot(freq3{i}, DGauss3{i})
    end
end
xlabel('Номер частоты')
ylabel('Полуширина, pix')

%%

%%