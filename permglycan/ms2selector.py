import pymzml
import pandas as pd
from bokeh.plotting import figure, output_file, show
from bokeh.models import Label, LabelSet, TeeHead, Arrow, Range1d
import pickle
from collections import OrderedDict
import os
import numpy as np
from threading import Thread
import glob

def threaded_print(*args):

    m = []
    m.extend(list(args))
    n = ' '.join(map(str, m))

    t = Thread(target=update_status, args=[n])
    t.start()
    t.join()


def update_status(s):
    print(s)

class thermo:
    def __init__(self, filename_in):

        self.filename = filename_in

        self.msrun = pymzml.run.Reader(filename_in, MSn_Precision=0.0005)


    def searchfragment(self, fragment_mz, rt_start, rt_end):

        s = os.path.splitext(self.filename)[0] + '_668frags.xlsx'
        writer = pd.ExcelWriter(s)

        # (scan number, scan time, prec mz, 668 mz, 668 intensity)
        result = []

        for index, spectrum in enumerate(self.msrun):
            # threaded_print('a. scan index = {}, level = {} at time = {}'.format(index,spectrum.ms_level, spectrum.scan_time))
            if spectrum.ms_level == 1: continue
            if spectrum.scan_time_in_minutes() < (rt_start): continue
            if spectrum.scan_time_in_minutes() > (rt_end): break

            matchList = spectrum.has_peak(fragment_mz)

            if matchList != []:
                # threaded_print('    b. matchlist = {}\n        c. RT = {} minutes'.format(matchList, spectrum.scan_time_in_minutes()))
                result.append([spectrum.ID,
                               spectrum.scan_time_in_minutes,
                               spectrum.selected_precursors[0]['mz'],
                               spectrum.selected_precursors[0]['charge'],
                               matchList[0][0],
                               matchList[0][1],
                               ])

        self.ms2hits = pd.DataFrame(result, columns=['scan', 'RT(minutes)', 'precursor_mz',
                                                     'precursor_charge', 'fragment_mz',
                                                     'fragment_intensity']
                                    )
        # threaded_print(self.ms2hits)



        self.ms2hits.to_excel(writer)
        writer.save()

def get_mzmls():
    z = glob.glob(r'O:\Analytics\Lab data backup\Orbitrap Fusion\2019\June\*.mzml')
    return z


if __name__ == '__main__':
    # # f = r'O:\Analytics\Lab data backup\Orbitrap Fusion\2019\June\perm glycans\matt_test.mzml'
    # f=r'O:\Analytics\Lab data backup\Orbitrap Fusion\2019\June\1N_NGlycan_Permethylated_IT_HCD25.mzML'
    # t = thermo(f)
    # t.searchfragment(668.3464, 1000/60,1200/60)
    list_of_files = get_mzmls()

    for f in list_of_files:
        threaded_print('starting file\n{}'.format(f))
        t = thermo(f)
        t.searchfragment(668.3464, 8, 49)