import numpy as np
from scipy.ndimage.filters import gaussian_filter
from astropy.io import fits
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.animation as manimation


def add_gaussian(ad_h, ad_dif, ad_crpix, ad_max_gauss, ad_a_gauss, ad_d_gauss):
    # Добавление h-ой гауссианы, полное количество равняется g на текущей
    # итерации.

    if np.any(np.isnan(ad_dif)):
        raise Exception(
            "Cannot fit " + str(ad_h) +
            "th gaussian according to this algorithm:"
            "ad_dif == NaN")

    # Поиск максимального значения в текущем массиве разности скана и
    # спокойного Солнца свернутого с ДНА + h-1 вписанных гауссиан с параметрами
    # max_gauss, a_gauss, d_gauss
    ad_max_gauss[ad_h] = np.argmax(ad_dif)  # координата максимума новой гауссианы

    # Выделение области в массиве dif, в которую будет вписана h-ая
    # гауссиана. Область заключает в себе найденный максимум max_gauss[h] и
    # ограничена справа и слева точками, в которых вторая производная меняет
    # знак.

    # Границы ищутся начиная с найденного максимума.
    right_bound = int(ad_max_gauss[ad_h])  # правая граница
    left_bound = int(ad_max_gauss[ad_h])  # левая граница

    # Поиск правой границы. На каждом шаге вправо проверяется знак второй
    # разностной производной.
    while (ad_dif[right_bound + 2] - ad_dif[right_bound + 1]) < \
            (ad_dif[right_bound + 1] - ad_dif[right_bound]):
        right_bound += 1

        # В случае, когда указатель выходит за правую границу массива,
        # выбрасывается исключение.
        if right_bound + 2 > len(ad_dif) - 1:
            raise Exception(
                "Cannot fit " + str(ad_h) +
                "th gaussian according to this algorithm:"
                "поиск приостановлен на правом краю массива")

    # Поиск левой границы
    while (ad_dif[left_bound - 2] - ad_dif[left_bound - 1]) < \
            (ad_dif[left_bound - 1] - ad_dif[left_bound]):
        left_bound -= 1

        # В случае, когда указатель выходит за левую границу массива,
        # выбрасывается исключение.
        if right_bound - 2 < 0:
            raise Exception(
                "Cannot fit " + str(ad_h) +
                "th gaussian according to this algorithm:"
                "поиск приостановлен на левом краю массива")

    # В том случае, если правая или левая граница не сдвинулась с начальной
    # позиции, кидает исключение.
    if right_bound == ad_max_gauss[ad_h] and left_bound == ad_max_gauss[ad_h]:
        raise Exception(
            "Cannot fit " + str(ad_h) +
            "th gaussian according to this algorithm:"
            "указатель не может быть сдвинут с начальной позиции")

    # Массив индексов области в массиве dif, в которую вписывается гауссиана.
    x_gauss = range(
        left_bound - int(ad_max_gauss[ad_h]),
        right_bound - int(ad_max_gauss[ad_h])
    )

    # Вычисление параметров гауссианы методом наименьших квадратов.
    a = np.zeros((len(x_gauss), 2))
    a[:, 0] += 1
    a[:, 1] += np.power(x_gauss, 2)
    b = np.log(
        np.abs(
            ad_dif[left_bound:right_bound]
        )
    )
    coef = np.dot(
        np.linalg.inv(
            np.dot(
                np.transpose(a), a
            )
        ),
        np.dot(
            np.transpose(a), b
        )
    )

    # Запись параметров гауссианы. Полуширина может оказаться чисто мнимой,
    # поэтому берется модуль.
    ad_a_gauss[ad_h] = np.exp(coef[0])  # амплитуда
    ad_d_gauss[ad_h] = np.abs(np.sqrt(-1 / (2 * coef[1])))

    x = np.array(range(int(-ad_crpix), len(ad_dif) - int(ad_crpix)))
    ad_dif -= ad_a_gauss[ad_h] * np.exp(
        - (x - ad_max_gauss[ad_h] + ad_crpix) ** 2 / (2 * ad_d_gauss[ad_h] ** 2)
    )

    if np.any(np.isnan(ad_d_gauss)):
        raise Exception(
            "Cannot fit " + str(ad_h) +
            "th gaussian according to this algorithm:"
            "ad_d_gauss == NaN")

    return ad_dif, ad_max_gauss, ad_a_gauss, ad_d_gauss


def minimize_quiet_sun(mi_crpix, mi_data, mi_convolution_template, mi_h,
                       mi_max_gauss, mi_a_gauss, mi_d_gauss):
    # Эта функция двигает шаблон спокойного солнца по временной координате,
    # чтобы постепенно уточнить время кульминции.

    # Сдвинутый шаблон спокойного Солнца
    mi_convolution = np.zeros(len(mi_data))
    # Середина сдвинутого шаблона спокойнго Солнца
    middle = np.round(len(mi_convolution) / 2).astype(int)
    # Середина шаблона спокойного Солнца
    middle_template = np.round(len(mi_convolution_template) / 2).astype(int)
    # Сдвиг
    shift = np.round(middle - mi_crpix).astype(int)

    mi_convolution = np.zeros(np.shape(mi_data))
    mi_convolution = np.interp(
        range(
            middle_template - middle + shift,
            middle_template + middle + shift
        ),
        range(
            len(mi_convolution_template)
        ),
        mi_convolution_template
    )

    # Сумма h гауссиан, которые были найдены на данной итерации.
    mi_gauss = np.zeros(len(mi_data))
    x = np.array(range(int(-mi_crpix), len(mi_data) - int(mi_crpix)))
    for f in range(mi_h):
        mi_gauss += mi_a_gauss[f] * np.exp(-(x - mi_max_gauss[f] + mi_crpix) ** 2 /
                                     (2 * mi_d_gauss[f] ** 2))

    # Подгонка шаблона спокойного Солнца к данным минус h гассуиан методом
    # наименьших квадратов.
    a = mi_convolution
    b = mi_data - mi_gauss
    mi_coef = 1 / np.dot(a, a) * np.dot(a, b)

    # Нахождение массива разности между сканом и сверткой спокойного Солнца
    # плюс h гассуиан.
    mi_dif = mi_data - mi_coef * mi_convolution - mi_gauss

    # Суммарная ошибка
    mi_discr = np.sum(np.abs(mi_dif))

    return mi_dif, mi_gauss, mi_discr, mi_coef * mi_convolution


def find_new_solar_center(fi_crpix, fi_data, fi_convolution_template, fi_h,
                          fi_max_gauss, fi_a_gauss, fi_d_gauss):
    # Нахождение нового солнечного центра методом наискорейшего спуска.
    # Меняется только центр Солнца crpix.

    new_crpix = fi_crpix  # Указатель на центр Солнца
    # Суммарная ошибка при текущем центре crpix
    _, _, cur_discr, _ = minimize_quiet_sun(
        fi_crpix, fi_data, fi_convolution_template, fi_h,
        fi_max_gauss, fi_a_gauss, fi_d_gauss
    )
    # Суммарная ошибка при сдвиге в положительном напривлении crpix + 1
    _, _, pos_step_discr, _ = minimize_quiet_sun(
        fi_crpix + 1, fi_data, fi_convolution_template, fi_h,
        fi_max_gauss, fi_a_gauss, fi_d_gauss
    )
    # Суммарная ошибка при сдвиге в отрицательном напривлении crpix - 1
    _, _, neg_step_discr, _ = minimize_quiet_sun(
        fi_crpix - 1, fi_data, fi_convolution_template, fi_h,
        fi_max_gauss, fi_a_gauss, fi_d_gauss
    )

    # Движение происходит только в одном изначально выбранном навравлении,
    # в сторону уменьшения суммарной ошибки.
    if pos_step_discr < cur_discr:
        while pos_step_discr < cur_discr:
            new_crpix += 1
            cur_discr = pos_step_discr
            _, _, pos_step_discr, _ = minimize_quiet_sun(
                new_crpix + 1, fi_data, fi_convolution_template, fi_h,
                fi_max_gauss, fi_a_gauss, fi_d_gauss
            )
    else:
        while neg_step_discr < cur_discr:
            new_crpix -= 1
            cur_discr = neg_step_discr
            _, _, neg_step_discr, _ = minimize_quiet_sun(
                new_crpix - 1, fi_data, fi_convolution_template, fi_h,
                fi_max_gauss, fi_a_gauss, fi_d_gauss
            )

    return new_crpix


def get_gaussian_parameters(
        get_data, get_crpix_in, get_num_gauss, get_convolution_template):
    #

    # Инициализация массивов, содержащих параметры вписанных гауссиан.
    get_max_gauss = np.zeros(get_num_gauss).astype(int)  # координаты максимумов
    # гауссиан
    get_a_gauss = np.zeros(get_num_gauss).astype(int)  # амплитуды гауссиан
    get_d_gauss = np.zeros(get_num_gauss).astype(int)  # полуширины гауссиан
    get_discr = 0  # Суммарная ошибка
    # Прочие буферные переменные.
    get_dif = np.zeros(len(get_data))  # массив остатков

    # Сглаживание пиков. Без этого работать не будет!
    get_data = gaussian_filter(get_data, 3)
    get_crpix = get_crpix_in

    for main_iter in range(get_num_gauss + 1):
        # На нулевом шаге (g = 0) очередная гауссиана не вписывается, функции
        # add_gaussian и minimize_quiet_sun не выполняются
        #crpix = crpix_in

        for get_h in range(main_iter):
            # В случае аварийной остановки функция возвращает массивы с
            # параметрами гауссиан, которые были записаны до возникновения
            # ошибки.
            try:
                # Вписывание h-ой гауссианы. Параметры новой гауссианы будут
                # записаны в массивы max_gauss, a_gauss, d_gauss на h позицию.
                get_dif, get_max_gauss, get_a_gauss, get_d_gauss = add_gaussian(
                    get_h, get_dif, get_crpix,
                    get_max_gauss, get_a_gauss, get_d_gauss
                )
            except Exception as e:
                print("Error:", e)

                return get_crpix, get_max_gauss, get_a_gauss, get_d_gauss, get_discr

            # Масштабирование шаблона спокойного Солнца по амплитуде после
            # добавление новой гауссианы.
            get_dif, get_gauss, get_discr, _ = minimize_quiet_sun(
                get_crpix, get_data, get_convolution_template, main_iter,
                get_max_gauss, get_a_gauss, get_d_gauss
            )

        # Изменение центра Солнца, эта переменная используется на
        # последующих шагах.
        #get_crpix = find_new_solar_center(
        #    get_crpix, get_data, get_convolution_template, main_iter,
        #    get_max_gauss, get_a_gauss, get_d_gauss)

        get_dif, get_gauss, get_discr, _ = minimize_quiet_sun(
            get_crpix, get_data, get_convolution_template, main_iter, get_max_gauss,
            get_a_gauss, get_d_gauss)

    return get_crpix, get_max_gauss, get_a_gauss, get_d_gauss, get_discr


def specify_gaussian_parameters(
    sp_data, sp_crpix_in, sp_num_gauss, sp_convolution_template,
    sp_max_gauss, sp_a_gauss, sp_d_gauss
):
#

    sp_discr = 0  # Суммарная ошибка
    # Прочие буферные переменные.
    sp_dif = np.zeros(len(sp_data))  # массив остатков

    # Сглаживание пиков. Без этого работать не будет
    sp_data = gaussian_filter(sp_data, 3)
    sp_crpix = sp_crpix_in

    #sp_crpix = find_new_solar_center(
    #    sp_crpix, sp_data, sp_convolution_template, sp_num_gauss,
    #    sp_max_gauss, sp_a_gauss, sp_d_gauss
    #)

    for i in range(sp_d_gauss.shape[0]):
        if np.isnan(sp_d_gauss[i]):
            sp_d_gauss[i] = 0
    _, sp_gauss, sp_discr, sp_quiet_sun = minimize_quiet_sun(
        sp_crpix, sp_data, sp_convolution_template, sp_num_gauss,
        sp_max_gauss, sp_a_gauss, sp_d_gauss
    )

    sp_dif = sp_data - sp_quiet_sun
    sp_max_gauss = np.zeros(sp_num_gauss)
    sp_a_gauss = np.zeros(sp_num_gauss)
    sp_d_gauss = np.zeros(sp_num_gauss)

    for sp_h in range(sp_num_gauss):
        # В случае аварийной остановки функция возвращает массивы с
        # параметрами гауссиан, которые были записаны до возникновения
        # ошибки.
        try:
            # Вписывание h-ой гауссианы. Параметры новой гауссианы будут
            # записаны в массивы max_gauss, a_gauss, d_gauss на h позицию.
            sp_dif, sp_max_gauss, sp_a_gauss, sp_d_gauss = add_gaussian(
                sp_h, sp_dif, sp_crpix, sp_max_gauss, sp_a_gauss, sp_d_gauss)
        except Exception as e:
            print("Error:", e)

            return sp_crpix, sp_max_gauss, sp_a_gauss, sp_d_gauss, sp_discr

    return sp_crpix, sp_max_gauss, sp_a_gauss, sp_d_gauss, sp_discr


def gaussian_sort(freq, max_gauss, a_gauss, d_gauss):
    # Сортировка гауссиан по положению максимума. Максимально допустимое
    # значение отклонения положения максимума между соседними частотами,
    # если отклонение меньше допустимого не найдено, то значение записывается
    # новый источник.

    new_max_gauss = np.zeros(np.shape(max_gauss))
    new_a_gauss = np.zeros(np.shape(a_gauss))
    new_d_gauss = np.zeros(np.shape(d_gauss))
    trig = 5
    num_gauss = np.shape(max_gauss)[1]
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


def get_spectr(fits_image_filename):
    hdul = fits.open(fits_image_filename)

    data = hdul[0].data
    CDELT1 = hdul[0].header['CDELT1']  # Коэффициент перевода из pix в arcsec
    R = int(round(hdul[0].header['SOLAR_R'] / CDELT1))  # Радиус Солнца (pix)
    CRPIX = hdul[0].header['CRPIX1']  # Коэффициент перевода из pix в arcsec

    freq = np.zeros(len(hdul[1].data))
    for i in range(len(freq)):
        freq[i] = hdul[1].data[i][3]

    convolution_template = pd.read_excel('quiet_sun_data.xlsx')
    num_gauss = 15
    max_gauss = np.zeros((len(freq), num_gauss))
    a_gauss = np.zeros((len(freq), num_gauss))
    d_gauss = np.zeros((len(freq), num_gauss))
    crpix = np.zeros(len(freq))
    discr = np.zeros(len(freq))

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

        for iter in range(10):
            crpix[freq_num],\
            max_gauss[freq_num, :],\
            a_gauss[freq_num, :], \
            d_gauss[freq_num, :],\
            discr[freq_num] = \
            specify_gaussian_parameters(
                data[freq_num, 0, :],
                crpix[freq_num],
                num_gauss,
                convolution_template_pix,
                max_gauss[freq_num, :],
                a_gauss[freq_num, :],
                d_gauss[freq_num, :]
            )
            #if freq_num == 17:
            #    one_plot(freq_num, data, convolution_template_pix, num_gauss,
            #             max_gauss,
            #             a_gauss, d_gauss, CDELT1, crpix, R, freq)

    make_video(crpix, CDELT1, data, convolution_template_pix, num_gauss, freq,
               max_gauss, a_gauss, d_gauss, R)

    max_gauss_sort, a_gauss_sort, d_gauss_sort, freq_sort = \
        gaussian_sort(freq, max_gauss, a_gauss, d_gauss)

    return max_gauss_sort, a_gauss_sort, d_gauss_sort, freq_sort, freq, crpix, \
           discr, R, CRPIX


def one_plot(freq_num, data, convolution_template_pix, num_gauss, max_gauss,
             a_gauss, d_gauss, CDELT1, crpix, R, freq):
    dif, gauss, discr, quiet_sun = minimize_quiet_sun(
        crpix[freq_num], data[freq_num, 0, :], convolution_template_pix,
        num_gauss,
        max_gauss[freq_num, :], a_gauss[freq_num, :], d_gauss[freq_num, :]
    )

    plot_x = np.linspace(
        - crpix[freq_num] * CDELT1,
        (data.shape[2] - crpix[freq_num]) * CDELT1,
        num=data.shape[2]
    )
    min_data = 0
    max_data = np.max(data[:, 0, :])

    plt.plot([0, 0], [min_data, max_data], 'k--')
    plt.plot([R * CDELT1, R * CDELT1], [min_data, max_data], 'k--')
    plt.plot([- R * CDELT1, - R * CDELT1], [min_data, max_data], 'k--')

    plt.plot(plot_x, data[freq_num, 0, :])
    plt.plot(plot_x, data[freq_num, 0, :] - dif)
    plt.plot(plot_x, quiet_sun)
    for j in range(num_gauss):
        one_gauss = a_gauss[freq_num, j] * np.exp(
            -(plot_x - (max_gauss[freq_num, j] - crpix[freq_num]) * CDELT1)
             ** 2 /
            (2 * (d_gauss[freq_num, j] * CDELT1) ** 2))
        plt.plot(plot_x, one_gauss)

    plt.axis([- 1.3 * R * CDELT1, 1.3 * R * CDELT1, 0, max_data])
    plt.xlabel('Distance from solar center, arcsec')
    plt.ylabel('Antenna temperature, K')
    plt.title('Gauss analisys')
    plt.text(- R * CDELT1, 0.9 * max_data, r'$\nu=$' + str(freq[freq_num]) +
             ' GHz')
    plt.show()


def make_video(crpix, CDELT1, data, convolution_template, h, freq,
               max_gauss, a_gauss, d_gauss, R):
    #

    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Movie Test', artist='Matplotlib',
                    comment='Movie support!')
    writer = FFMpegWriter(fps=2, bitrate=-1, metadata=metadata)

    fig = plt.figure()

    with writer.saving(fig, "writer_test.mp4", 100):
        for i in range(len(freq)):
            dif, gauss, discr, quiet_sun = minimize_quiet_sun(
                crpix[i], data[i, 0, :], convolution_template, h,
                max_gauss[i, :], a_gauss[i, :], d_gauss[i, :]
            )

            plot_x = np.linspace(
                - crpix[i] * CDELT1,
                (data.shape[2] - crpix[i]) * CDELT1,
                num=data.shape[2]
            )
            min_data = - 0.1 * np.max(data[:, 0, :])
            max_data = np.max(data[:, 0, :])

            plt.plot([0, 0], [min_data, max_data], 'k--')
            plt.plot([R * CDELT1, R * CDELT1], [min_data, max_data], 'k--')
            plt.plot([- R * CDELT1, - R * CDELT1], [min_data, max_data], 'k--')

            plt.plot(plot_x, data[i, 0, :])
            plt.plot(plot_x, data[i, 0, :] - dif)
            plt.plot(plot_x, quiet_sun)
            plt.plot(plot_x, dif, 'k-')
            for j in range(h):
                one_gauss = a_gauss[i, j] * np.exp(
                    -(plot_x - (max_gauss[i, j] - crpix[i]) * CDELT1)
                     ** 2 /
                    (2 * (d_gauss[i, j] * CDELT1) ** 2))
                plt.plot(plot_x, one_gauss)

            plt.axis([- 1.3 * R * CDELT1, 1.3 * R * CDELT1, min_data, max_data])
            plt.xlabel('Distance from solar center, arcsec')
            plt.ylabel('Antenna temperature, K')
            plt.title('Gauss analisys')
            plt.text(- R * CDELT1, 0.9 * max_data, r'$\nu=$' + str(freq[i]) +
                     ' GHz')

            writer.grab_frame()
            fig.clf()


fits_image_filename = '20150715_121933_sun0_out.fits'  # Имя файла

max_gauss_sort, a_gauss_sort, d_gauss_sort, freq_sort, freq, crpix, discr, R, CRPIX\
    = \
    get_spectr(fits_image_filename)

for i in range(len(freq_sort)):
    plt.plot([freq[0], freq[-1]], [CRPIX, CRPIX], 'k--')
    plt.plot([freq[0], freq[-1]], [CRPIX + R, CRPIX + R], 'k--')
    plt.plot([freq[0], freq[-1]], [CRPIX - R, CRPIX - R], 'k--')
    plt.plot(freq_sort[i], max_gauss_sort[i])
plt.xlabel('Frequency, GHz')
plt.ylabel('Position, pix')
plt.title('Position')
plt.show()

for i in range(len(freq_sort)):
    plt.plot(freq_sort[i], a_gauss_sort[i])
plt.xlabel('Frequency, GHz')
plt.ylabel('Amplitude, K')
plt.title('Amplitude')
plt.show()

for i in range(len(freq_sort)):
    plt.plot(freq_sort[i], d_gauss_sort[i])
plt.xlabel('Frequency, GHz')
plt.ylabel('FWHM, pix')
plt.title('Full width at half maximum, pix')
plt.show()

plt.plot(freq, discr)
plt.xlabel('Frequency, GHz')
plt.ylabel('Error')
plt.title('Total error')
plt.show()

plt.plot(freq, crpix)
plt.show()
