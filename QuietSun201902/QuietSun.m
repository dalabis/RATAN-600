clear
close all

% Усредненное спокойное Солнце будет состоять из отцентрованных сканов
% длиной 4096 пикселей
quietSun(1:4096, 1:84) = 0;

fileNames = {'20190115_122255_sun+0_out_edit.fits', ...
             '20190225_122643_sun+0_out_edit.fits', ...
             '20190227_122623_sun+0_out_edit.fits', ...
             '20190228_122613_sun+0_out_edit.fits', ...
             '20190301_122602_sun+0_out_edit.fits', ...
             '20190302_122550_sun+0_out_edit.fits', ...
             '20190303_122538_sun+0_out_edit.fits'};
for fileName = fileNames
    data = fitsread(char(fileName));
    % Нужно обрезать моменты включения генератора шума
    cutDataLeft = 200;
    cutDataRight = 200;
    data = data(:,1+cutDataLeft:end-cutDataRight,:);
    % Чтение значения положения солнечного из заголовка fits. файла
    info = fitsinfo(char(fileName));
    CRPIX1Success = 0;
    i = 1;
    while (i <= size(info.PrimaryData.Keywords, 1)) && ~CRPIX1Success
        if isequal(info.PrimaryData.Keywords{i,1}, 'CRPIX1')
            CRPIX1 = info.PrimaryData.Keywords{i,2};
            CRPIX1Success = 1;
        end
        
        i = i + 1;
    end
    
    for i = 1:size(data, 3)
        center = fix(2048 + (size(data,2)/2 - (CRPIX1 - cutDataLeft)));
        
        if mod(size(data,2),2) == 0
            quietSun(center-fix(size(data,2)/2):center+fix(size(data,2)/2)-1,i) = ...
                quietSun(center-fix(size(data,2)/2):center+fix(size(data,2)/2)-1,i) + averaging(data(1,:,i))';
        else
            quietSun(center-fix(size(data,2)/2):center+fix(size(data,2)/2),i) = ...
                quietSun(center-fix(size(data,2)/2):center+fix(size(data,2)/2),i) + averaging(data(1,:,i))';
        end
    end
end

figure
hold on
axis tight

for i = 1:size(data, 3)
    plot(quietSun(:,i) ./ size(fileNames, 2))
end

% Запись шаблона в fits. файл
fitswrite(quietSun ./ size(fileNames, 2),'template.fits');

%% Averaging
% Эта функция строит непараметрическую регрессию с нормальным ядром
% Так пики получаются более гладкими
% Но ширину ядра нужно подбирать для каждой частоты отдельно (Можно связать с шириной ДНА, но почему то не работает)
function averData = averaging(data)
    % Ширина окна
    h = 10;
    % Число отчетов
    m = length(data);
    %
    x = 1:m;
    
    % Усреднение
    averData(1:length(h),1:m) = 0;
    for i = 1:length(h)
        for j = 1:m
            K = 1 / (h(i) * sqrt(2 * pi)) .* exp( - ( j - x ).^2 ./ ( 2 * h(i)^2 ) );
            averData(i, j) = sum( K .* data ) / sum(K);
        end
    end
end