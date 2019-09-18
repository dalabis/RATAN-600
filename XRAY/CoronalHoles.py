import numpy as np
from PIL import Image
from astropy.io import fits
import cv2 as cv


def low_intensity_regions_map(hdul):
    # Алгоритм считает гистограмму в каждом подизображении размера 1/2, 1/3,
    # 1/4, 1/5. На каждой гистограмме ищется локальный минимум с интенсивностью
    # 0.3 - 0.7 от максимального значения QS_intensity (Quiet Sun intensity).
    # Ширина минимума должна быть не менее 6 столбцов гистограммы. Далее
    # алгоритм строит гистогрмму минимумов и находит максимум. Области с
    # интенсивность меньшей чем найденая по гистограмме отмечаются как области
    # с пониженной интенсивностью LIRs (Low Intensity Regions).

    x0 = int(round(hdul[0].header['CRPIX1']))  # Координаты центра Солнца (pix)
    y0 = int(round(hdul[0].header['CRPIX2']))
    CDELT1 = hdul[0].header['CDELT1']  # Коэффициент перевода из pix в arcsec
    R = int(round(hdul[0].header['RSUN'] / CDELT1))  # Радиус Солнца (pix)

    dataSphere = hdul[0].data

    # Исходные декартовы координаты пикселей на диске Солнца. За пределами диска
    # координаты равны 0
    x = np.zeros((2 * R, 2 * R)).astype(int)
    y = np.zeros((2 * R, 2 * R)).astype(int)
    for i in range(2 * R):
        for j in range(2 * R):
            if (i - R) ** 2 + (j - R) ** 2 < R ** 2:
                x[i, j] = j - R
                y[i, j] = i - R

    # Исходные сферические координаты пикселей на поверхности Солнца. Широта и
    # долгота отсчитываются от ближайшей точки на поверхности Солнца к наблюдателю
    # (центр диска). За пределами диска координаты равны 0.
    phi = np.arcsin(y / R)  # Широта
    lamb = np.arcsin(x / (R * np.cos(phi)))  # Долгота

    # Переход к новым декартовым координатам. (Долгота, Широта) сферические
    # -> (Долгота, Широта) декартовы.
    X = R * lamb  # Долгота
    Y = R * np.sin(phi)  # Широта

    # Перевод долготы и широты в индексы нового массива
    X = np.round((X - np.amin(X)) / (np.amax(X) - np.amin(X)) * 2 * R - 1).astype(
        int)
    Y = np.round(Y - np.amin(Y)).astype(int)

    # Формирование нового массива изображения в декартовых координатах
    dataCart = np.zeros((2 * R, 2 * R))

    # Заполнение нового массива. Некоторые пиксели останутся 0. Далее алгоритм 0
    # пиксели будет игнорировать
    for i in range(2 * R):
        for j in range(2 * R):
            dataCart[X[i, j], Y[i, j]] = dataSphere[x0 - x[i, j], y0 - y[i, j]]

    poolMin = np.zeros(
        3 ** 2 + 5 ** 2 + 7 ** 2 + 9 ** 2)  # Массив найденных минимумов
    l = 0  # Счетчик минимумов
    k = 0  # Счетчик для поиска локального минимума на гистограмме
    boolMin = 1  # Проверка на локальный минимум
    dataMax = np.round(np.amax(dataCart)).astype(int)  # Максимальная интенсивность
    value = 0  # Итенсивность текущего элемента изображения

    # Гитограмма интенсивностей всего изображения. Пиксели с 0 интенсивностью не
    # учитываются
    hist_QS_intensity, _ = np.histogram(dataCart, range(1, dataMax))

    # Максимальная интенсивность
    QS_intensity = np.argmax(hist_QS_intensity)

    for subNum in range(2, 6):
        subLen = np.round(2 * R / subNum).astype(int)  # Длина края подизображения

        # Изображения берутся 'внахлест' с шагом subLen/2
        for i in range(2 * (subNum - 1)):
            for j in range(2 * (subNum - 1)):
                poolImage = dataCart[np.round(subLen * i / 2).astype(int):
                                     subLen + np.round(subLen * i / 2).astype(int),
                                     np.round(subLen * j / 2).astype(int):
                                     subLen + np.round(subLen * j / 2).astype(
                                int)]  # Подизображение

                # Построение гистограммы интенсивностей локального подизображения.
                # Пиксели с 0 интенсивностью не учитываются.
                hist_pool_intensity, _ = np.histogram(poolImage, range(1, dataMax))

                # Минимум ищется в пределах 0.3 - 0.7 QS_intensity
                k = np.round(0.3 * QS_intensity).astype(int)

                while k < 0.7 * QS_intensity:
                    # Поиск локального минимума, т.е. элемента, который является
                    # минимальным в своей окресности. В данном случае окрестность
                    # из шести столбцов гистограммы

                    boolMin = 1

                    for p in range(-5, 5):
                        if hist_pool_intensity[k] > hist_pool_intensity[k + p]:
                            boolMin = 0
                            break

                    if boolMin and (poolMin[l] == 0 or poolMin[l] > hist_pool_intensity[k]):
                        poolMin[l] = k

                    k += 1

                l += 1

    # нулевые значения игнорируются
    hist_trash, _ = np.histogram(poolMin, range(1, dataMax))

    # Пиксели со значением интенсивности ниже LIR_thrash являютя областями с
    # пониженной интенсивностью LIRs.
    LIR_thrash = np.argmax(hist_trash)

    # Бинарный массив, значение 1 говорит о наличии области пониженной
    # интенсивности.
    LIR_binary = np.zeros(np.shape(dataSphere))
    for i in range(2 * R):
        for j in range(2 * R):
            if dataCart[X[i, j], Y[i, j]] < LIR_thrash and \
                    dataCart[X[i, j], Y[i, j]] != LIR_thrash:
                LIR_binary[x0 - x[i, j], y0 - y[i, j]] = 1

    LIR_binary = LIR_binary.astype('uint8')
    contours, _ = cv.findContours(LIR_binary,
                                  cv.RETR_TREE,
                                  cv.CHAIN_APPROX_SIMPLE)

    return contours


fits_image_filename = 'saia_00193_fd_20190705_232152.fts'  # Имя файла
hdul = fits.open(fits_image_filename)

low_intensity_regions = low_intensity_regions_map(hdul)

# Формирование цветного RGB изображения
image = Image.fromarray(hdul[0].data)
image = image.convert('RGB')
image = np.asarray(image, dtype="int32")

# Добавление границ корональных дыр на изображение
image = cv.drawContours(cv.UMat(image),
                        low_intensity_regions,
                        -1,
                        (0, 0, 255),
                        1)

# Сохранение изображения
cv.imwrite("test.jpg", image)

hdul.close()