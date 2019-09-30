import numpy as np
from scipy.ndimage.filters import gaussian_filter
from astropy.io import fits
import pandas as pd
from matplotlib import pyplot as plt


def add_gaussian(h, dif, max_gauss, a_gauss, d_gauss):
    # Добавление h-ой гауссианы, полное количество равняется g на текущей
    # итерации.

    # Поиск максимального значения в текущем массиве разности скана и
    # спокойного Солнца свернутого с ДНА + h-1 вписанных гауссиан с параметрами
    # max_gauss, a_gauss, d_gauss
    max_gauss[h] = np.argmax(dif)  # координата максимума новой гауссианы

    # Выделение области в массиве dif, в которую будет вписана h-ая
    # гауссиана. Область заключает в себе найденный максимум max_gauss[h] и
    # ограничена справа и слева точками, в которых вторая производная меняет
    # знак.

    # Границы ищутся начиная с найденного максимума.
    right_bound = max_gauss[h]  # правая граница
    left_bound = max_gauss[h]  # левая граница

    # Поиск правой границы. На каждом шаге вправо проверяется знак второй
    # разностной производной.
    while (dif[right_bound + 2] - dif[right_bound + 1]) < \
            (dif[right_bound + 1] - dif[right_bound]):
        right_bound += 1

        # В случае, когда указатель выходит за правую границу массива,
        # выбрасывается исключение.
        if right_bound + 2 > len(dif) - 1:
            raise Exception("Cannot fit " + h + "th gaussian according to "
                                                "this algorithm: поиск "
                                                "приостановлен на правом краю "
                                                "массива")

    # Поиск левой границы
    while (dif[left_bound - 2] - dif[left_bound - 1]) < \
            (dif[left_bound - 1] - dif[left_bound]):
        left_bound -= 1

        # В случае, когда указатель выходит за левую границу массива,
        # выбрасывается исключение.
        if right_bound - 2 < 0:
            raise Exception("Cannot fit " + h + "th gaussian according to "
                                                "this algorithm: поиск "
                                                "приостановлен на левом краю "
                                                "массива")

    # В том случае, если правая или левая граница не сдвинулась с начальной
    # позиции, кидает исключение.
    if right_bound == max_gauss[h] or left_bound == max_gauss[h]:
        raise Exception("Cannot fit " + h + "th gaussian according to "
                                            "this algorithm: указатель не "
                                            "может быть сдвинут с начальной "
                                            "позиции")

    # Массив индексов области в массиве dif, в которую вписывается гауссиана.
    x_gauss = range(left_bound - max_gauss[h], right_bound - max_gauss[h])

    # Вычисление параметров гауссианы методом наименьших квадратов.
    a = np.zeros((len(x_gauss), 2))
    a[:, 0] += 1
    a[:, 1] += np.power(x_gauss, 2)
    b = np.log(np.abs(dif[left_bound:right_bound]))
    coef = np.dot(np.linalg.inv(np.dot(np.transpose(a), a)),
                  np.dot(np.transpose(a), b))

    # Запись параметров гауссианы. Полуширина может оказаться чисто мнимой,
    # поэтому берется модуль.
    a_gauss[h] = np.exp(coef[0])  # амплитуда
    d_gauss[h] = np.abs(np.sqrt(-1 / (2 * coef[1])))

    return max_gauss, a_gauss, d_gauss


def minimize_quiet_sun(crpix, data, convolution_template, h,
                       max_gauss, a_gauss, d_gauss):
    # Эта функция двигает шаблон спокойного солнца по временной координате,
    # чтобы постепенно уточнить время кульминции.

    # Сдвинутый шаблон спокойного Солнца
    convolution = np.zeros(len(data))
    # Середина сдвинутого шаблона спокойнго Солнца
    middle = np.round(len(convolution) / 2).astype(int)
    # Середина шаблона спокойного Солнца
    middle_template = np.round(len(convolution_template) / 2).astype(int)
    # Сдвиг
    shift = np.round(middle - crpix).astype(int)

    convolution = np.interp(range(middle_template - middle + shift,
                                  middle_template + middle + shift),
                            range(len(convolution_template)),
                            convolution_template)

    # Сумма h гауссиан, которые были найдены на данной итерации.
    gauss = np.zeros(len(data))
    x = np.array(range(int(-crpix), len(data) - int(crpix)))
    for f in range(h):
        gauss += a_gauss[f] * np.exp(-(x - max_gauss[f] + crpix) ** 2 /
                                     (2 * d_gauss[f] ** 2))

    # Подгонка шаблона спокойного Солнца к данным минус h гассуиан методом
    # наименьших квадратов.
    a = convolution
    b = data - gauss
    coef = 1 / np.dot(a, a) * np.dot(a, b)

    # Нахождение массива разности между сканом и сверткой спокойного Солнца
    # плюс h гассуиан.
    dif = data - coef * convolution - gauss

    # Суммарная ошибка
    discr = np.sum(np.abs(dif))

    return dif, gauss, discr


def find_new_solar_center(crpix, data, convolution_template, h,
                          max_gauss, a_gauss, d_gauss):
    # Нахождение нового солнечного центра методом наискорейшего спуска.
    # Меняется только центр Солнца crpix.

    new_crpix = crpix  # Указатель на центр Солнца
    # Суммарная ошибка при текущем центре crpix
    _, _, cur_discr = minimize_quiet_sun(crpix, data, convolution_template, h,
                                         max_gauss, a_gauss, d_gauss)
    # Суммарная ошибка при сдвиге в положительном напривлении crpix + 1
    _, _, pos_step_discr = minimize_quiet_sun(crpix + 1, data,
                                              convolution_template, h,
                                              max_gauss, a_gauss, d_gauss)
    # Суммарная ошибка при сдвиге в отрицательном напривлении crpix - 1
    _, _, neg_step_discr = minimize_quiet_sun(crpix - 1, data,
                                              convolution_template, h,
                                              max_gauss, a_gauss, d_gauss)

    # Движение происходит только в одном изначально выбранном навравлении,
    # в сторону уменьшения суммарной ошибки.
    if pos_step_discr < cur_discr:
        while pos_step_discr < cur_discr:
            new_crpix += 1
            cur_discr = pos_step_discr
            _, _, pos_step_discr = minimize_quiet_sun(new_crpix + 1, data,
                                                      convolution_template,
                                                      h, max_gauss, a_gauss,
                                                      d_gauss)
    else:
        while neg_step_discr < cur_discr:
            new_crpix -= 1
            cur_discr = neg_step_discr
            _, _, neg_step_discr = minimize_quiet_sun(new_crpix - 1, data,
                                                      convolution_template,
                                                      h, max_gauss, a_gauss,
                                                      d_gauss)

    return new_crpix


def get_gaussian_parameters(data, crpix_in, num_gauss, convolution_template):
    #

    # Инициализация массивов, содержащих параметры вписанных гауссиан.
    max_gauss = np.zeros(num_gauss).astype(int)  # координаты максимумов
    # гауссиан
    a_gauss = np.zeros(num_gauss).astype(int)  # амплитуды гауссиан
    d_gauss = np.zeros(num_gauss).astype(int)  # полуширины гауссиан
    discr = 0  # Суммарная ошибка
    # Прочие буферные переменные.
    dif = np.zeros(len(data))  # массив остатков

    # Сглаживание пиков. Без этого работать не будет!
    data = gaussian_filter(data, 5)

    for g in range(0, num_gauss + 1):
        # На нулевом шаге (g = 0) очередная гауссиана не вписывается, функции
        # add_gaussian и minimize_quiet_sun не выполняются
        crpix = crpix_in

        for h in range(0, g):
            # В случае аварийной остановки функция возвращает массивы с
            # параметрами гауссиан, которые были записаны до возникновения
            # ошибки.
            try:
                # Вписывание h-ой гауссианы. Параметры новой гауссианы будут
                # записаны в массивы max_gauss, a_gauss, d_gauss на h позицию.
                max_gauss, a_gauss, d_gauss = add_gaussian(h, dif, max_gauss,
                                                           a_gauss, d_gauss)
            except Exception as e:
                print("Error:", e)

                return crpix, max_gauss, a_gauss, d_gauss, discr

            # Масштабирование шаблона спокойного Солнца по амплитуде после
            # добавление новой гауссианы.
            dif, gauss, discr = minimize_quiet_sun(crpix, data,
                                                   convolution_template, g-1,
                                                   max_gauss, a_gauss, d_gauss)

        # Изменение центра Солнца, эта переменная используется на
        # последующих шагах.
        crpix = find_new_solar_center(crpix, data, convolution_template, g,
                                      max_gauss, a_gauss, d_gauss)

        dif, gauss, discr = minimize_quiet_sun(crpix, data,
                                               convolution_template, g,
                                               max_gauss, a_gauss, d_gauss)

    return crpix, max_gauss, a_gauss, d_gauss, discr


def gaussian_sort(freq, max_gauss, a_gauss, d_gauss):
    # Сортировка гауссиан по положению максимума. Максимально допустимое
    # значение отклонения положения максимума между соседними частотами,
    # если отклонение меньше допустимого не найдено, то значение записывается
    # новый источник.

    new_max_gauss = np.zeros(np.shape(max_gauss))
    new_a_gauss = np.zeros(np.shape(a_gauss))
    new_d_gauss = np.zeros(np.shape(d_gauss))
    trig = 10
    size = num_gauss
    new_size = num_gauss
    min = 0
    min_ind = 0

    new_max_gauss[0, :] = max_gauss[0, :]
    new_a_gauss[0, :] = a_gauss[0, :]
    new_d_gauss[0, :] = d_gauss[0, :]

    for i in range(1, len(freq)):
        size = new_size

        for j in range(num_gauss):
            min = trig
            min_ind = 0

            for l in range(i - 1):
                for k in range(size):
                    if np.abs(new_max_gauss[l, k] - max_gauss[i, j]) < min:
                        min = np.abs(new_max_gauss[l, k] - max_gauss[i, j])
                        min_ind = k

            if min == trig:
                # Если не удалось найти источник на предыдущей частоте,
                # который расположен не дальше чем trig пикселей, то для
                # источника добавляется новая строка.
                new_size += 1
                new_max_gauss = np.column_stack((new_max_gauss, np.transpose(
                    np.zeros(len(freq)))))
                new_a_gauss = np.column_stack((new_a_gauss, np.transpose(
                    np.zeros(len(freq)))))
                new_d_gauss = np.column_stack((new_d_gauss, np.transpose(
                    np.zeros(len(freq)))))
                new_max_gauss[i, new_size - 1] = max_gauss[i, j]
                new_a_gauss[i, new_size - 1] = a_gauss[i, j]
                new_d_gauss[i, new_size - 1] = d_gauss[i, j]
            else:
                new_max_gauss[i, min_ind] = max_gauss[i, j]
                new_a_gauss[i, min_ind] = a_gauss[i, j]
                new_d_gauss[i, min_ind] = d_gauss[i, j]

    list_max_gauss = []
    list_a_gauss = []
    list_d_gauss = []
    list_freq = []

    for i in range(np.shape(new_max_gauss)[1]):
        object_max_gauss = []
        object_a_gauss = []
        object_d_gauss = []
        object_freq = []

        for j in range(np.shape(new_max_gauss)[0]):
            if new_max_gauss[j, i] != 0:
                object_max_gauss.append(new_max_gauss[j, i])
                object_a_gauss.append(new_a_gauss[j, i])
                object_d_gauss.append(new_d_gauss[j, i])
                object_freq.append(freq[j])

        list_max_gauss.append(object_max_gauss)
        list_a_gauss.append(object_a_gauss)
        list_d_gauss.append(object_d_gauss)
        list_freq.append(object_freq)

    # Ограничение на колличество точек в спектре
    spectr_len = 15
    max_gauss_out = []
    a_gauss_out = []
    d_gauss_out = []
    freq_out = []

    for spectr in list_max_gauss:
        if len(spectr) > spectr_len:
            max_gauss_out.append(spectr)
    for spectr in list_a_gauss:
        if len(spectr) > spectr_len:
            a_gauss_out.append(spectr)
    for spectr in list_d_gauss:
        if len(spectr) > spectr_len:
            d_gauss_out.append(spectr)
    for spectr in list_freq:
        if len(spectr) > spectr_len:
            freq_out.append(spectr)

    return max_gauss_out, a_gauss_out, d_gauss_out, freq_out


fits_image_filename = '20150715_121933_sun0_out.fits'  # Имя файла
hdul = fits.open(fits_image_filename)

data = hdul[0].data
CDELT1 = hdul[0].header['CDELT1']  # Коэффициент перевода из pix в arcsec
R = int(round(hdul[0].header['SOLAR_R'] / CDELT1))  # Радиус Солнца (pix)
CRPIX = hdul[0].header['CRPIX1']  # Коэффициент перевода из pix в arcsec

freq = np.zeros(len(hdul[1].data))
for i in range(len(freq)):
    freq[i] = hdul[1].data[i][3]

convolution_template = pd.read_excel('quiet_sun_data.xlsx')
num_gauss = 8
max_gauss = np.zeros((len(freq), num_gauss))
a_gauss = np.zeros((len(freq), num_gauss))
d_gauss = np.zeros((len(freq), num_gauss))
crpix = np.zeros(len(freq))
discr = np.zeros(len(freq))

print(CRPIX)

for freq_num in range(len(freq)):
    freq_template = 1

    for i in range(1, len(convolution_template.columns)):
        if np.abs(float(convolution_template.columns[i]) - freq[freq_num]) <\
                np.abs(float(convolution_template.columns[freq_template]) -
                       freq[freq_num]):
            freq_template = i

    convolution_template_arcsec = np.zeros(2 * len(convolution_template) - 1)
    x_arcsec = np.zeros(2 * len(convolution_template) - 1)

    convolution_template_arcsec[len(convolution_template)-1:] = \
        np.array(convolution_template[
                     convolution_template.columns[freq_template]
                 ])
    convolution_template_arcsec[:len(convolution_template)] = \
        np.array(convolution_template[
                     convolution_template.columns[freq_template]
                 ])[::-1]
    x_arcsec[len(convolution_template)-1:] = \
        np.array(convolution_template[convolution_template.columns[0]])
    x_arcsec[:len(convolution_template)] = \
        - np.array(convolution_template[convolution_template.columns[0]])[::-1]
    # Перевод из угловых секунд в пиксели
    x_pix = (x_arcsec / 750 * (R * CDELT1)) / CDELT1
    # Интерполяция шаблона по пикселям
    x_pix_interp = np.array(range(-2000, 2000))
    convolution_template_pix = np.interp(x_pix_interp, x_pix,
                                         convolution_template_arcsec)

    crpix[freq_num] = CRPIX

    crpix[freq_num], max_gauss[freq_num, :], a_gauss[freq_num, :],\
        d_gauss[freq_num, :], discr[freq_num] =\
        get_gaussian_parameters(data[freq_num, 0, :], crpix[freq_num],
                                num_gauss, convolution_template_pix)

max_gauss_sort, a_gauss_sort, d_gauss_sort, freq_sort = \
    gaussian_sort(freq, max_gauss, a_gauss, d_gauss)

for i in range(len(freq_sort)):
    plt.plot(freq_sort[i], max_gauss_sort[i])
plt.show()

for i in range(len(freq_sort)):
    plt.plot(freq_sort[i], a_gauss_sort[i])
plt.show()

for i in range(len(freq_sort)):
    plt.plot(freq_sort[i], d_gauss_sort[i])
plt.show()
