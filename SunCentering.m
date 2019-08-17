function [CRPIX, maxGauss, AGauss, DGauss, currentDiscr] = SunCentering(data, CRPIX, R, ii, numGauss, alpha, Dfreq, sunShape)

    t = -R:R;
    % Инициальзация массивов, которые содержат параметры вписнных
    % гауссиан
    maxGauss(1:numGauss) = 0;
    AGauss(1:numGauss) = 0;
    DGauss(1:numGauss) = 0;
    
    for g = 0:numGauss
        % Инициальзация массивов, которые содержат параметры вписнных
        % гауссиан
        % На нулевом шаге (g = 0) гауссиана не вписывается, т.е. блоки Add
        % и Min не выполняются

        for h = 1:g  
            %% Add
            % Добавление h-ой гауссианы, полное колличество g на текущем
            % шаге
            
            % Поиск максимального значения в текущем массиве разности
            % скана и спокойного Солнца свернутого с ДНА + h-1 гауссиан с
            % текущими параметрами maxGauss, AGauss, DGauss
            [~, maxGauss(h)] = max(dif);
            
            % Сглаживание шума
            dif = smooth(dif,5)';
            
            % Выделение области в массиве dif для аппроксимации очередной
            % h-ой гауссианой
            % Область заключает в себе максимум и ограничена точками
            % перегиба первого порядка (хорошо работает)
            rightBound = maxGauss(h);
            leftBound = maxGauss(h);
            % На каждом шаге вправо проверяется знак второй разностной
            % производной
            % !!! Возможен еще один тип исключения: Subscript indices must 
            % either be real positive integers or logicals (нужно это
            % учесть)
            % Поиск правой границы
            while (dif(rightBound + 2) - dif(rightBound + 1)) < (dif(rightBound + 1) - dif(rightBound))
                rightBound = rightBound + 1;
            end
            % Поиск левой границы
            while (dif(leftBound - 2) - dif(leftBound - 1)) < (dif(leftBound - 1) - dif(leftBound))
                leftBound = leftBound - 1;
            end
            % В том случае, если правая или левая граница не сдвинулась,
            % кидает исклчение, которое нужно обработать в основной функции
            if rightBound == maxGauss(h) || leftBound == maxGauss(h)
                disp('Cannot fit th Gaussian according to this algorithm.')
                
                return
            end
            
            % Чтобы избежать отрицательной разницы, массив dif 
            % "приподнимается" так, чтобы гауссиана была вписана с
            % макимальной амплитудой
            rightMin = maxGauss(h);
            leftMin = maxGauss(h);
            % Поиск правой границы
            while dif(rightMin + 1) < dif(rightMin)
                rightMin = rightMin + 1;
            end
            % Поиск левой границы
            while dif(leftMin + 1) < dif(leftMin)
                leftMin = leftMin + 1;
            end
            % Определение минимального значения
            if dif(rightMin) < dif(leftMin)
                minGauss = dif(rightMin);
            else
                minGauss = dif(leftMin);
            end
            % Если минимальное значение больше 0, то сдвига не должно быть
            if minGauss > 0
                minGauss = 0;
            end
            
            % Индексы области в массиве dif, в которой вписывается
            % гауссиана
            xGauss = leftBound-maxGauss(h):rightBound-maxGauss(h);
            
            % Гауссиана вписывается методом наименьших квадратов
            A1(1:length(xGauss)) = 1;
            A2 = xGauss.^2;
            A = [A1', A2'];
            A1 = [];
            b = log(abs(dif(leftBound:rightBound) - minGauss));
            coef = ( A' * A ) \ ( A' * b' );
            % Параметры очередной гауссианы записываются в массив
            % параметров
            AGauss(h) = exp(coef(1));
            DGauss(h) = sqrt( -1 / ( 2 * coef(2) ) );

            %% Min
            % Подгонка свертки спокойного Солнца с ДНА + h текущих гауссиан
            % по двум параметрам:
            %   1. Сдвиг по оси времени
            %   2. Масштабирование по амплитуде (домножение на константу)
            % Методом наименьших квадратов с регулизационным параметром 
            % alpha (входной параметр)
            
            % Свертка спокойного Солнца (уже свернутого с вертикальной ДНА)
            % с горизонтальной диаграммой направленности
            % Горизонтальная диаграмма направленности антенны
            antennaPattern = @(x) exp( - x.^2 / ( 2 * (Dfreq(ii)/2.355)^2 ) );
            % Свертка
            convolution(1:length(x)) = 0;
            for j = 1:length(x)
                convolution(j) = trapz(sunShape(:,ii)' .* antennaPattern(x(j)-t));
            end
            
            % Сумма всех h текущих гауссиан
            gauss(1:size(data, 2)) = 0;
            for f = 1:h
                gauss = gauss + AGauss(f).*exp(-(x-maxGauss(f)+CRPIX).^2./(2.*DGauss(f).^2));
            end

            % Подгонка свертки Спокойного солнсца с ДНА к скану - сумма
            % всех текущих гауссиан
            % МНК с гребневой регуляризацией (alpha параметр регуляризации)
            A = convolution';
            b = data(1,:,ii) - gauss;
            coef = (A'*A + alpha) \ (A'*b');

            % нахождение разницы между сканом и сверткой спокойного Солнца
            % с ДНА + h гауссиан
            dif = data(1,:,ii) - coef.*convolution - gauss;
        end

        %% Find
        % Нахождение нового Солнечного центра методом наискорейшего спуска
        % На этом шаге меняется только центр покойного Солнца
        % Функция discrepancy вычисляет относиельную разницу между сверткой
        % спокойного Солнца ДНА + h гауссиан и одномерным сканом
        newCRPIX = CRPIX;
        currentDiscr = discrepancy(CRPIX, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
        posStepDiscr = discrepancy(CRPIX+1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
        negStepDiscr = discrepancy(CRPIX-1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
        
        % Движение происходит только в одном изначально выбранном
        % направлении
        % Движение вправо, пока относительная разница не перестает меняться
        if posStepDiscr < currentDiscr
            while posStepDiscr < currentDiscr
                newCRPIX = newCRPIX + 1;
                currentDiscr = posStepDiscr;
                posStepDiscr = discrepancy(newCRPIX+1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
            end
        % Движение влево, пока относительная разница не перестает меняться
        elseif negStepDiscr < currentDiscr
            while negStepDiscr < currentDiscr
                newCRPIX = newCRPIX - 1;
                currentDiscr = negStepDiscr;
                negStepDiscr = discrepancy(newCRPIX-1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
            end
        end
        
        % Изменение центра Солнца, эта переменная используется на
        % последующих шагах
        CRPIX = newCRPIX;

        %% Plot
        % plot quiet Sun shape and RATAN scan with first gaussian
        x = 1-CRPIX:1:size(data,2)-CRPIX;

        antennaPattern = @(x) exp( - x.^2 / ( 2 * (Dfreq(ii)/2.355)^2 ) );
        convolution(1:length(x)) = 0;

        for j = 1:length(x)
            convolution(j) = trapz(sunShape(:,ii)' .* antennaPattern(x(j)-t));
        end

        %figure
        %hold on
        
        gauss(1:size(data, 2)) = 0;
        for h = 1:g
            %plot(x, AGauss(h).*exp(-(x-maxGauss(h)+CRPIX).^2./(2.*DGauss(h).^2)))
            gauss = gauss + AGauss(h).*exp(-(x-maxGauss(h)+CRPIX).^2./(2.*DGauss(h).^2));
        end
        A = convolution';
        b = data(1,:,ii) - gauss;
        coef = (A'*A + alpha) \ (A'*b');
        %[~,point] = min(data(1,:,ii) - convolution - gauss);
        %coef = 1/convolution(point)*data(1,point,ii);

        % diferences between scan and quiet Sun shape
        dif = data(1,:,ii) - coef.*convolution;
        
        %plot(x, data(1,:,ii))
        %plot(x, coef.*convolution + gauss, '--')
        %plot(x, coef.*convolution)
        %plot(x, dif - gauss, '-')
    end
end

%% sun shape discrepancy function %%%%%%%%%%%%%%%%%%%%%%%%%%
% this function finds discrepancy from RATAN scan and quiet Sun shape
function discr = discrepancy(CRPIX, R, Dfreq, data, sunShape, ii, maxGauss, AGauss, DGauss, alpha)
    x = 1-CRPIX:1:size(data,2)-CRPIX;
    t = -R:R;
    antennaPattern = @(x) exp( - x.^2 / ( 2 * (Dfreq(ii)/2.355)^2 ) );
    convolution(1:length(x)) = 0;
    
    % layer with gaussians
    gauss(1:length(x)) = 0;
    for j = 1:length(maxGauss)
        gauss = gauss + AGauss(j).*exp(-(x-maxGauss(j)+CRPIX).^2./(2.*DGauss(j).^2));
    end
    
    for j = 1:length(x)
        convolution(j) = trapz(sunShape' .* antennaPattern(x(j)-t));
    end
    
    A = convolution';
    b = data(1,:,ii) - gauss;
    coef = (A'*A + alpha) \ (A'*b');
    %[~,point] = min(data(1,:,ii) - convolution - gauss);
    %coef = 1/convolution(point)*data(1,point,ii);
    
    dif = data(1,:,ii) - coef.*convolution - gauss;
    
    discr = sum(abs(dif));
end