import pymzml
import pandas as pd
from bokeh.plotting import figure, output_file, show
from bokeh.models import Label, LabelSet, TeeHead, Arrow, Range1d
import pickle
from collections import OrderedDict
import os
import numpy as np


class thermo:
    def __init__(self, filename_in):

        self.msrun = pymzml.run.Reader(filename_in)

        # need to create an index of runs so that the msrun can be sliced by RT.

        # self.rt_series = pd.DataFrame([[spectrum['scan time'],
        #                               spectrum['ms level'],
        #                               spectrum.highest_peaks(50)] \
        #                                for spectrum in self.msrun],
        #                               columns=['rt','level','peaks'])


    def xic(self, mz, rt, rt_window):
        # returns a time series of the summed intensities of all the mz's parsed
        # as a list of tuples [(rt, intensity)...]

        timeDependentIntensities = []
        matchList=[[0,0,0]]
        for index, spectrum in enumerate(self.msrun):
            print('a. scan index = {}'.format(index))
            if spectrum.ms_level != 1: continue
            if spectrum.scan_time[0] < (rt - rt_window): continue
            if spectrum.scan_time[0] > (rt + rt_window): break

            matchList = spectrum.has_peak(mz)

            if matchList != []:
                print('b. matchlist = {}\nc. RT = {} minutes'.format(matchList, spectrum.scan_time[0]))
                for mz, I in matchList:
                    timeDependentIntensities.append([spectrum.scan_time[0], I, mz])
            else:
                timeDependentIntensities.append([spectrum.scan_time[0], 0.0, mz])

        return pd.DataFrame(timeDependentIntensities, columns=['rt', 'intensity', 'mz'])


        # for rt, i, mz in timeDependentIntensities:
        #     print('{0:5.3f} {1:13.4f}       {2:10}'.format(rt, i, mz))


def fix_obo_in_mzml_file(file_in):
    with open(
            file_in,
            'rb')as f:
        a = f.read().replace(b'23:06:2017', b'4.0.1')

    with open(
            file_in,
            'wb') as f:
        f.write(a)

    return file_in


def get_list_of_ions(gradient='long'):
    if gradient=='short':
        result = [(674.3197, 9.5, 3, 'N1H1_3+', 'blue'),
                  (1156.5235, 9.5,3, 'N1H1S1_2+', 'red'),
                  (771.3515, 9.5, 3, 'N1H1S1_3+', 'darkorange'),
                  (1302.0712, 9.5, 3, 'N1H1S2_2+', 'lime'),
                  (868.3833, 9.5, 3, 'N1H1S2_3+', 'limegreen'),
                  (990.0940, 9.5, 3, 'N2H2S1_3+', 'black'),
                  (519.5791, 12, 2, 'AVK-PAK_2+', 'gray'),
                  ]
    elif gradient=='long':
        result = [(552.6089, 19, 3, 'naked_3+', 'purple'),
                  (674.3197, 22, 5, 'N1H1_3+', 'blue'),
                  (1156.5235, 22, 5, 'N1H1S1_2+', 'red'),
                  (771.3515, 22, 5, 'N1H1S1_3+', 'darkorange'),
                  (1302.0712, 22, 5, 'N1H1S2_2+', 'lime'),
                  (868.3833, 26, 3, 'N1H1S2_3+', 'limegreen'),
                  (990.0940, 26, 3, 'N2H2S1_3+', 'black'),
                  # (560.2589, 57, 2, 'VYF-FDR_2+', 'gray'),
                  ]

    return result


def do_xic(filename_in, pickle_save_name):
    """

    :param filename_in:
    :param pickle_save_name:
    :return:
    """
    # filename_in = r'O:\Analytics\Lab data backup\Orbitrap Fusion\2019\Apr\Calnexin\sCANX_New_Ctrl_Try+GluC_24ug_4ug_ul_C18_OT_EtHCD40_targeted_mz200_2000_90min_RD7_2.mzML'

    # this replaces an obsolete OBO entry with a valid one; needed for MSConverted mzml files.
    filename_in = fix_obo_in_mzml_file(filename_in)

    combined_dfs = OrderedDict()

    # initialise the plot
    p = figure(plot_width=400, plot_height=400)

    # cycle through a list of ions obtained from a function and get the XIC
    # this is an inefficient way as we have to reopen the MS as you cycle through each ion
    # would be better to open once and grab all the xic needed.
    for ion_tuple in get_list_of_ions():
        t = thermo(filename_in=filename_in)

        # get the XIC based on mz, rt, rt_window
        a = t.xic(ion_tuple[0], ion_tuple[1], ion_tuple[2])
        combined_dfs[str(ion_tuple[0])] = a
    #     p.line(a['rt'], a['intensity'], line_width=2, legend=ion_tuple[3], color=ion_tuple[4])
    #
    # p.legend.location = "top_right"
    # p.legend.click_policy = "hide"
    # output_file('test.html')
    # show(p)

    with open(pickle_save_name, 'wb') as f:
        pickle.dump(combined_dfs, f)

    return combined_dfs


def get_default_plot(title='Extracted Ion Chromatograms of V[42]EDSKPDTTAPPSSPK[57] O-glycopeptides'):
    p = figure(plot_width=800, plot_height=400, title=title or 'Plot')
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    # p.y_range = Range1d(-2.5e6, 2.5e6)
    return p


def plot_xic_from_pickle(pkl_file_name):

    with open(pkl_file_name, 'rb') as f:
        combined_dict_of_dfs = pickle.load(f)

    p = get_default_plot()

    for key_name, ion_tuple in zip(combined_dict_of_dfs, get_list_of_ions()):
        a = combined_dict_of_dfs[key_name]
        p.line(a['rt'], a['intensity'], line_width=2, legend=ion_tuple[3], color=ion_tuple[4])


    output_file(os.path.splitext(pkl_file_name)[0] + '.html')
    p.legend.location = "top_right"
    p.legend.click_policy = "hide"
    show(p)


def plot_butterfly_xic_from_pickles(pkl_file_name1, pkl_file_name2, rt_correction=(0,0)):

    with open(pkl_file_name1, 'rb') as f:
        combined_dict_of_dfs1 = pickle.load(f)

        p = get_default_plot()

        AUC1 = 0

        for key_name, ion_tuple in zip(combined_dict_of_dfs1, get_list_of_ions()):
            a = combined_dict_of_dfs1[key_name]
            p.line(a['rt']+rt_correction[0], a['intensity'], line_width=2, legend=ion_tuple[3], color=ion_tuple[4],
                   )

            AUC1 += np.trapz(y=a['intensity'])

    with open(pkl_file_name2, 'rb') as f:
        combined_dict_of_dfs2 = pickle.load(f)
        AUC2 = 0
        for key_name, ion_tuple in zip(combined_dict_of_dfs2, get_list_of_ions()):
            a = combined_dict_of_dfs2[key_name]
            p.line(a['rt']+rt_correction[1], -a['intensity'], line_width=2, legend=ion_tuple[3], color=ion_tuple[4])
            AUC2 += np.trapz(y=a['intensity'])

    output_file('butterfly.html')
    p.output_backend = 'canvas'
    p.legend.location = "top_right"
    p.legend.click_policy = "hide"
    p.yaxis.axis_label = "Intensity"
    p.xaxis.axis_label = "Retention Time (minutes)"

    # add annotation and labels
    p.add_layout(Label(x=6.5, y=-2e6, x_units='screen',
                         text='DOX Induced', render_mode='css')
                 )
    p.add_layout(Label(x=6.5, y=2e6, x_units='screen',
                       text='Control', render_mode='css')
                 )
    p.add_layout(Label(x=700, y=-2e6, x_units='screen',
                       text='AUC(Control/DOX)={:.4f}'.format(AUC1/AUC2), render_mode='css',
                       text_align='right')
                 )
    # p.add_layout(Label(x=7.9, y=1.4e6,
    #                    text='0-NeuAc', render_mode='css',
    #                    text_align='center',
    #                    text_color='blue')
    #              )
    # p.add_layout(Label(x=8.7, y=1.4e6,
    #                    text='1-NeuAc', render_mode='css',
    #                    text_align='center',
    #                    text_color='orange')
    #              )
    # p.add_layout(Label(x=9.9, y=1.4e6,
    #                    text='2-NeuAc & Extended', render_mode='css',
    #                   text_align='left',
    #                    text_color='green')
    #              )
    #

    # p.add_layout(Arrow(start=TeeHead(line_color='blue'), end=TeeHead(line_color='blue'),
    #                    x_start=7.8, x_end=8.25,
    #                    y_start=1e6, y_end=1e6,
    #                    line_width=1, line_color="blue",
    #                    ))
    # p.add_layout(Arrow(start=TeeHead(line_color='orange'), end=TeeHead(line_color='orange'),
    #                    x_start=8.4, x_end=8.85,
    #                    y_start=1e6, y_end=1e6,
    #                    line_color='orange'
    #                    ))
    # p.add_layout(Arrow(start=TeeHead(line_color='green'), end=TeeHead(line_color='green'),
    #                    x_start=9.3, x_end=10.2,
    #                    y_start=1e6, y_end=1e6,
    #                    line_color='green'
    #                    ))



    show(p)

def execute(only_plot_from_pickle=True):

    ctrl = 'control_xic_data_802.pkl'
    dox = 'dox_xic_data_902.pkl'
    DIRPATH = r'O:\Analytics\Lab data backup\Orbitrap Fusion\2019\Apr\CalnexinLongRunMay'

    if not only_plot_from_pickle:

        combined_dfs_control = do_xic(
            os.path.join(DIRPATH, r'802_sCANX_Ctrl_TryGluC_12ug_C18_ITHCD30_mz200_1650_210min_RD7_2.mzML'),
            ctrl)

        combined_dfs_dox = do_xic(
            os.path.join(DIRPATH, r'902_sCANX_DOX_TryGluC_12ug_C18_ITHCD30_mz200_1650_210min_RD8_2.mzML'),
            dox)

    # plot_xic_from_pickle('control_xic_data.pkl')
    # plot_xic_from_pickle('dox_xic_data.pkl')
    plot_butterfly_xic_from_pickles(ctrl, dox, rt_correction=(0,-0.5))

if __name__ == "__main__":

    execute(only_plot_from_pickle=True)