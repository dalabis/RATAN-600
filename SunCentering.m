function [CRPIX, maxGauss, AGauss, DGauss, currentDiscr] = ...
    SunCentering(data, CRPIX, freqNum, numGauss, alpha, convolutionTemplate, graphType, v)

    % ������������� �������� ���������� ��������� ��������� ��������
    maxGauss(1:numGauss) = 0;
    AGauss(1:numGauss) = 0;
    DGauss(1:numGauss) = 0;
    
    % ����������� �����
    data(1,:,freqNum) = averaging(data(1,:,freqNum));
    data(1,:,freqNum) = data(1,:,freqNum) - min(data(1,:,freqNum));
    
    for g = 0:numGauss
        % ������������� ��������, ������� �������� ��������� ��������
        % ��������
        % �� ������� ���� (g = 0) ��������� �� �����������, �.�. ����� Add
        % � Min �� �����������

        for h = 1:g  
            %% Add
            % ���������� h-�� ���������, ������ ����������� g �� �������
            % ��������
            
            % ����� ������������� �������� � ������� ������� ��������
            % ����� � ���������� ������ ���������� � ��� + h-1 �������� �
            % �������� ����������� maxGauss, AGauss, DGauss
            [~, maxGauss(h)] = max(dif);
            
            
            % ��������� ������� � ������� dif ��� ������������� ���������
            % h-�� ����������
            % ������� ��������� � ���� �������� � ���������� �������
            % �������� ������� ������� (������ ��������)
            rightBound = maxGauss(h);
            leftBound = maxGauss(h);
            % �� ������ ���� ������ ����������� ���� ������ ����������
            % �����������
            % !!! �������� ��� ���� ��� ����������: Subscript indices must 
            % either be real positive integers or logicals (����� ���
            % ������)
            % ����� ������ �������
            while (dif(rightBound + 2) - dif(rightBound + 1)) < (dif(rightBound + 1) - dif(rightBound))
                rightBound = rightBound + 1;
                
                if rightBound + 2 > length(dif)
                    disp(['Cannot fit ', num2str(h), 'th Gaussian according to this algorithm: ����� ������������� �� ������ ���� �������'])
                    
                    return
                end
            end
            % ����� ����� �������
            while (dif(leftBound - 2) - dif(leftBound - 1)) < (dif(leftBound - 1) - dif(leftBound))
                leftBound = leftBound - 1;
                
                if leftBound - 2 < 1
                    disp(['Cannot fit ', num2str(h), 'th Gaussian according to this algorithm: ����� ������������� �� ����� ���� �������'])
                    
                    return
                end
            end
            % � ��� ������, ���� ������ ��� ����� ������� �� ����������,
            % ������ ���������, ������� ����� ���������� � �������� �������
            if rightBound == maxGauss(h) || leftBound == maxGauss(h)
                disp(['Cannot fit ', num2str(h), 'th Gaussian according to this algorithm: ������������� ���������'])
                
                return
            end
            
            % ������� ������� � ������� dif, � ������� �����������
            % ���������
            xGauss = leftBound-maxGauss(h):rightBound-maxGauss(h);
            
            % ��������� ����������� ������� ���������� ���������
            A1(1:length(xGauss)) = 1;
            A2 = xGauss.^2;
            A = [A1', A2'];
            A1 = [];
            b = log(abs(dif(leftBound:rightBound)));
            coef = ( A' * A ) \ ( A' * b' );
            % ��������� ��������� ��������� ������������ � ������
            % ����������
            AGauss(h) = exp(coef(1));
            % ���������� ����� ��������� ����� ������, ������� �������
            % ������
            DGauss(h) = abs(sqrt( -1 / ( 2 * coef(2) ) ));

            %% Min
            % �������� ������� ���������� ������ � ��� + h ������� ��������
            % �� ���� ����������:
            %   1. ����� �� ��� �������
            %   2. ��������������� �� ��������� (���������� �� ���������)
            % ������� ���������� ��������� � ��������������� ���������� 
            % alpha (������� ��������)
            
            % ������ ���������� ������
            % ������� ��� � ��������� ��������� ������� ��� ������
            % ���������� ������
            x = 1-CRPIX:1:size(data,2)-CRPIX;
            convolution(1:length(x)) = 0;
            middle = round(length(convolution)/2);
            middleTemplate = round(length(convolutionTemplate)/2);
            shift = round(middle - CRPIX);
            
            if mod(length(x),2) == 1 && mod(length(middleTemplate),2) == 1
                convolution = convolutionTemplate(middleTemplate-middle+1+shift:middleTemplate+middle-1+shift);
            else
                convolution = convolutionTemplate(middleTemplate-middle+1+shift:middleTemplate+middle+shift);
            end
            
            % ����� ���� h ������� ��������
            gauss(1:size(data, 2)) = 0;
            for f = 1:h
                gauss = gauss + AGauss(f).*exp(-(x-maxGauss(f)+CRPIX).^2./(2.*DGauss(f).^2));
            end

            % �������� ������� ���������� ������� � ��� � ����� - �����
            % ���� ������� ��������
            % ��� � ��������� �������������� (alpha �������� �������������)
            A = convolution';
            b = data(1,:,freqNum) - gauss;
            coef = (A'*A + alpha) \ (A'*b');

            % ���������� ������� ����� ������ � �������� ���������� ������
            % � ��� + h ��������
            dif = data(1,:,freqNum) - coef.*convolution - gauss;
        end

        %% Find
        % ���������� ������ ���������� ������ ������� ������������� ������
        % �� ���� ���� �������� ������ ����� ��������� ������
        % ������� discrepancy ��������� ������������ ������� ����� ��������
        % ���������� ������ ��� + h �������� � ���������� ������
        newCRPIX = CRPIX;
        currentDiscr = discrepancy(CRPIX, data, convolutionTemplate, freqNum, maxGauss, AGauss, DGauss, alpha);
        posStepDiscr = discrepancy(CRPIX+1, data, convolutionTemplate, freqNum, maxGauss, AGauss, DGauss, alpha);
        negStepDiscr = discrepancy(CRPIX-1, data, convolutionTemplate, freqNum, maxGauss, AGauss, DGauss, alpha);
        
        % �������� ���������� ������ � ����� ���������� ���������
        % �����������
        % �������� ������, ���� ������������� ������� �� ��������� ��������
        if posStepDiscr < currentDiscr
            while posStepDiscr < currentDiscr
                newCRPIX = newCRPIX + 1;
                currentDiscr = posStepDiscr;
                posStepDiscr = discrepancy(newCRPIX+1, data, convolutionTemplate, freqNum, maxGauss, AGauss, DGauss, alpha);
            end
        % �������� �����, ���� ������������� ������� �� ��������� ��������
        elseif negStepDiscr < currentDiscr
            while negStepDiscr < currentDiscr
                newCRPIX = newCRPIX - 1;
                currentDiscr = negStepDiscr;
                negStepDiscr = discrepancy(newCRPIX-1, data, convolutionTemplate, freqNum, maxGauss, AGauss, DGauss, alpha);
            end
        end
        
        % ��������� ������ ������, ��� ���������� ������������ ��
        % ����������� �����
        CRPIX = newCRPIX;

        %% Plot
        % plot quiet Sun shape and RATAN scan with first gaussian
        % ������ ���������� ������
        % ������� ��� � ��������� ��������� ������� ��� ������
        % ���������� ������
        x = 1-CRPIX:1:size(data,2)-CRPIX;
        convolution(1:length(x)) = 0;
        middle = round(length(convolution)/2);
        middleTemplate = round(length(convolutionTemplate)/2);
        shift = round(middle - CRPIX);
        
        if mod(length(x),2) == 1 && mod(length(middleTemplate),2) == 1
            convolution = convolutionTemplate(middleTemplate-middle+1+shift:middleTemplate+middle-1+shift);
        else
            convolution = convolutionTemplate(middleTemplate-middle+1+shift:middleTemplate+middle+shift);
        end
        
        if isequal(graphType, 'AllSteps')
            figure
            hold on

            gauss(1:size(data, 2)) = 0;
            for h = 1:g
                plot(x, AGauss(h).*exp(-(x-maxGauss(h)+CRPIX).^2./(2.*DGauss(h).^2)))
                gauss = gauss + AGauss(h).*exp(-(x-maxGauss(h)+CRPIX).^2./(2.*DGauss(h).^2));
            end
            A = convolution';
            b = data(1,:,freqNum) - gauss;
            coef = (A'*A + alpha) \ (A'*b');

            % diferences between scan and quiet Sun shape
            dif = data(1,:,freqNum) - coef.*convolution;

            plot(x, data(1,:,freqNum))
            plot(x, coef.*convolution + gauss, '--')
            plot(x, coef.*convolution)
            plot(x, dif - gauss, '-')
        else
            gauss(1:size(data, 2)) = 0;
            for h = 1:g
                gauss = gauss + AGauss(h).*exp(-(x-maxGauss(h)+CRPIX).^2./(2.*DGauss(h).^2));
            end
            A = convolution';
            b = data(1,:,freqNum) - gauss;
            coef = (A'*A + alpha) \ (A'*b');

            % diferences between scan and quiet Sun shape
            dif = data(1,:,freqNum) - coef.*convolution;
            
            if h == numGauss
                figure
                hold on
                
                for h = 1:g
                    plot(x, AGauss(h).*exp(-(x-maxGauss(h)+CRPIX).^2./(2.*DGauss(h).^2)))
                end
                
                plot(x, data(1,:,freqNum))
                plot(x, coef.*convolution + gauss, '--')
                plot(x, coef.*convolution)
                plot(x, dif - gauss, '-')
                
                xlabel('Distance from Solar center, arcsec')
                ylabel('Antenna temperature, K')
                
                axis tight
                
                % �������� �����
                frame = getframe(gcf);
                
                % ������ ����� � ����������
                for j = 1:10
                    writeVideo(v,frame);
                end
                
                % �������� ������������ ����
                close gcf
            end
        end
    end
end

%% Discrepancy
% ��� ������� ������� ������������� ������� ����� �������� ����������
% ������ � ��� + h �������� � ���������� ������
function discr = discrepancy(CRPIX, data, convolutionTemplate, freqNum, maxGauss, AGauss, DGauss, alpha)
    % ������ ���������� ������
    % ������� ��� � ��������� ��������� ������� ��� ������
    % ���������� ������
    x = 1-CRPIX:1:size(data,2)-CRPIX;
    convolution(1:length(x)) = 0;
    middle = round(length(convolution)/2);
    middleTemplate = round(length(convolutionTemplate)/2);
    shift = round(middle - CRPIX);
    
    if mod(length(x),2) == 1 && mod(length(middleTemplate),2) == 1
        convolution = convolutionTemplate(middleTemplate-middle+1+shift:middleTemplate+middle-1+shift);
    else
        convolution = convolutionTemplate(middleTemplate-middle+1+shift:middleTemplate+middle+shift);
    end
    
    % layer with gaussians
    gauss(1:length(x)) = 0;
    for j = 1:length(maxGauss)
        gauss = gauss + AGauss(j).*exp(-(x-maxGauss(j)+CRPIX).^2./(2.*DGauss(j).^2));
    end
    
    A = convolution';
    b = data(1,:,freqNum) - gauss;
    coef = (A'*A + alpha) \ (A'*b');
    
    dif = data(1,:,freqNum) - coef.*convolution - gauss;
    
    discr = sum(abs(dif));
end

%% Averaging
% ��� ������� ������ ����������������� ��������� � ���������� �����
% ��� ���� ���������� ����� ��������
% �� ������ ���� ����� ��������� ��� ������ ������� �������� (����� ������� � ������� ���)
function averData = averaging(data)
    % ������ ����
    h = 5;
    % ����� �������
    m = length(data);
    %
    x = 1:m;
    
    % ����������
    averData(1:length(h),1:m) = 0;
    for i = 1:length(h)
        for j = 1:m
            K = 1 / (h(i) * sqrt(2 * pi)) .* exp( - ( j - x ).^2 ./ ( 2 * h(i)^2 ) );
            averData(i, j) = sum( K .* data ) / sum(K);
        end
    end
end