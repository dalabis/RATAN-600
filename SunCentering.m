function [CRPIX, maxGauss, AGauss, DGauss, currentDiscr] = SunCentering(data, CRPIX, R, ii, numGauss, alpha, Dfreq, sunShape)

    t = -R:R;
    % ������������� ��������, ������� �������� ��������� ��������
    % ��������
    maxGauss(1:numGauss) = 0;
    AGauss(1:numGauss) = 0;
    DGauss(1:numGauss) = 0;
    
    for g = 0:numGauss
        % ������������� ��������, ������� �������� ��������� ��������
        % ��������
        % �� ������� ���� (g = 0) ��������� �� �����������, �.�. ����� Add
        % � Min �� �����������

        for h = 1:g  
            %% Add
            % ���������� h-�� ���������, ������ ����������� g �� �������
            % ����
            
            % ����� ������������� �������� � ������� ������� ��������
            % ����� � ���������� ������ ���������� � ��� + h-1 �������� �
            % �������� ����������� maxGauss, AGauss, DGauss
            [~, maxGauss(h)] = max(dif);
            
            % ����������� ����
            dif = smooth(dif,5)';
            
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
            end
            % ����� ����� �������
            while (dif(leftBound - 2) - dif(leftBound - 1)) < (dif(leftBound - 1) - dif(leftBound))
                leftBound = leftBound - 1;
            end
            % � ��� ������, ���� ������ ��� ����� ������� �� ����������,
            % ������ ���������, ������� ����� ���������� � �������� �������
            if rightBound == maxGauss(h) || leftBound == maxGauss(h)
                disp('Cannot fit th Gaussian according to this algorithm.')
                
                return
            end
            
            % ����� �������� ������������� �������, ������ dif 
            % "��������������" ���, ����� ��������� ���� ������� �
            % ����������� ����������
            rightMin = maxGauss(h);
            leftMin = maxGauss(h);
            % ����� ������ �������
            while dif(rightMin + 1) < dif(rightMin)
                rightMin = rightMin + 1;
            end
            % ����� ����� �������
            while dif(leftMin + 1) < dif(leftMin)
                leftMin = leftMin + 1;
            end
            % ����������� ������������ ��������
            if dif(rightMin) < dif(leftMin)
                minGauss = dif(rightMin);
            else
                minGauss = dif(leftMin);
            end
            % ���� ����������� �������� ������ 0, �� ������ �� ������ ����
            if minGauss > 0
                minGauss = 0;
            end
            
            % ������� ������� � ������� dif, � ������� �����������
            % ���������
            xGauss = leftBound-maxGauss(h):rightBound-maxGauss(h);
            
            % ��������� ����������� ������� ���������� ���������
            A1(1:length(xGauss)) = 1;
            A2 = xGauss.^2;
            A = [A1', A2'];
            A1 = [];
            b = log(abs(dif(leftBound:rightBound) - minGauss));
            coef = ( A' * A ) \ ( A' * b' );
            % ��������� ��������� ��������� ������������ � ������
            % ����������
            AGauss(h) = exp(coef(1));
            DGauss(h) = sqrt( -1 / ( 2 * coef(2) ) );

            %% Min
            % �������� ������� ���������� ������ � ��� + h ������� ��������
            % �� ���� ����������:
            %   1. ����� �� ��� �������
            %   2. ��������������� �� ��������� (���������� �� ���������)
            % ������� ���������� ��������� � ��������������� ���������� 
            % alpha (������� ��������)
            
            % ������� ���������� ������ (��� ���������� � ������������ ���)
            % � �������������� ���������� ��������������
            % �������������� ��������� �������������� �������
            antennaPattern = @(x) exp( - x.^2 / ( 2 * (Dfreq(ii)/2.355)^2 ) );
            % �������
            convolution(1:length(x)) = 0;
            for j = 1:length(x)
                convolution(j) = trapz(sunShape(:,ii)' .* antennaPattern(x(j)-t));
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
            b = data(1,:,ii) - gauss;
            coef = (A'*A + alpha) \ (A'*b');

            % ���������� ������� ����� ������ � �������� ���������� ������
            % � ��� + h ��������
            dif = data(1,:,ii) - coef.*convolution - gauss;
        end

        %% Find
        % ���������� ������ ���������� ������ ������� ������������� ������
        % �� ���� ���� �������� ������ ����� ��������� ������
        % ������� discrepancy ��������� ������������ ������� ����� ��������
        % ���������� ������ ��� + h �������� � ���������� ������
        newCRPIX = CRPIX;
        currentDiscr = discrepancy(CRPIX, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
        posStepDiscr = discrepancy(CRPIX+1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
        negStepDiscr = discrepancy(CRPIX-1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
        
        % �������� ���������� ������ � ����� ���������� ���������
        % �����������
        % �������� ������, ���� ������������� ������� �� ��������� ��������
        if posStepDiscr < currentDiscr
            while posStepDiscr < currentDiscr
                newCRPIX = newCRPIX + 1;
                currentDiscr = posStepDiscr;
                posStepDiscr = discrepancy(newCRPIX+1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
            end
        % �������� �����, ���� ������������� ������� �� ��������� ��������
        elseif negStepDiscr < currentDiscr
            while negStepDiscr < currentDiscr
                newCRPIX = newCRPIX - 1;
                currentDiscr = negStepDiscr;
                negStepDiscr = discrepancy(newCRPIX-1, R, Dfreq, data, sunShape(:,ii), ii, maxGauss, AGauss, DGauss, alpha);
            end
        end
        
        % ��������� ������ ������, ��� ���������� ������������ ��
        % ����������� �����
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