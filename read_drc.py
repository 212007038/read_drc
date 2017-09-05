###############################################################################
# Python script used to read in Datex-Ohmeda Record file and parse out
# ECG12 sub-records.

###############################################################################

# region Import region
import os
import argparse  # command line parser
from struct import *
import numpy as np

# endregion

# region Globals region

__version__ = '1.0'  # version of script

WF_FORMAT = '<hBBBB hBBBB hBBBB hBBBB hBBBB hBBBB hBBBB hBBBB HHHH'

# Record main-type
DRC_HEADER_SIZE = 40
DRI_MT_PHDB = 0
DRI_MT_WAVE = 1
DRI_MT_ARMS = 2
DRI_MT_B2B_CTRL = 3
DRI_MT_ALARM = 4
DRI_MT_NETWORK = 5
DRI_MT_PRINT = 6
DRI_MT_SERV = 7
DRI_MT_FO = 8
DRI_MT_ACK = 9
DRI_MT_EVENT = 10
DRI_MT_MEM = 11
DRI_MT_MODES = 12
DRI_MT_IND = 13
DRI_MT_TEST = 14
DRI_MT_INPH = 15
DRI_MT_ADU_INTF_PHDB = 16
DRI_MT_DEBUG = 17
DRI_MT_SNAPSHOT = 18
DRI_MT_CASE = 19


# Subrecord ID Value Subrecord Chapter DRI level
DRI_WF_CMD = 0  # Waveform request
DRI_WF_ECG1 = 1  # Old ECG waveforms
DRI_WF_ECG2 = 2  # Old ECG waveforms
DRI_WF_ECG3 = 3  # Old ECG waveforms
DRI_WF_INVP1 = 4  # Invasive Blood Pressure Waveforms
DRI_WF_INVP2 = 5  # Invasive Blood Pressure Waveforms
DRI_WF_INVP3 = 6  # Invasive Blood Pressure Waveforms
DRI_WF_INVP4 = 7  # Invasive Blood Pressure Waveforms
DRI_WF_PLETH = 8  # Pleth Waveform
DRI_WF_CO2 = 9  # CO2 Waveform
DRI_WF_O2 = 10  # O2 Waveform
DRI_WF_N2O = 11  # N2O Waveform
DRI_WF_AA = 12  # Anesthesia Agent Waveform
DRI_WF_AWP = 13  # PAW Waveform
DRI_WF_FLOW = 14  # Flow Waveform
DRI_WF_RESP = 15  # Respiration Waveform
DRI_WF_INVP5 = 16  # Invasive Blood Pressure Waveforms DRI_LEVEL_97->
DRI_WF_INVP6 = 17  # Invasive Blood Pressure Waveforms DRI_LEVEL_97->
DRI_WF_EEG1 = 18  # EEG Waveform DRI_LEVEL_98->
DRI_WF_EEG2 = 19  # EEG Waveform DRI_LEVEL_98->
DRI_WF_EEG3 = 20  # EEG Waveform DRI_LEVEL_98->
DRI_WF_EEG4 = 21  # EEG Waveform DRI_LEVEL_98->
DRI_WF_ECG12 = 22  # ECG 12-lead Waveform DRI_LEVEL_98->
DRI_WF_VOL = 23  # Airway volume Waveform DRI_LEVEL_98->
DRI_WF_TONO_PRESS = 24  # Tonometry pressure Waveform DRI_LEVEL_98->
DRI_WF_VENT_AWP = 25  # It is never used. DRI_LEVEL_98->
DRI_WF_VENT_FLOW = 26  # It is never used. DRI_LEVEL_98->
DRI_WF_VENT_VOL = 27  # It is never used. DRI_LEVEL_98->
DRI_WF_NIBP = 28  # NIBP cuff pressure waveform. DRI_LEVEL_99->
DRI_WF_SPI_LOOP_STATUS = 29  # Spirometry loop DRI_LEVEL_99->


datex_hdr_fmt = '<HBBHLBBHH24s'
ecg12_wf_data_fmt = '<hBBBB hBBBB hBBBB hBBBB hBBBB hBBBB hBBBB hBBBB HHHH'
wf_hdr_fmt = '<HHH HHHHHHH'
suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']

# endregion


def humansize(nbytes):
    i = 0
    while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
    f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
    return '%s %s' % (f, suffixes[i])


def build_sign_array():
    # Build a gigantic multipication array
    m_array = np.empty([65536, 16], dtype=np.int16)
    l = [0x8000, 0x2000, 0x0800, 0x0200, 0x0080, 0x0020, 0x0008, 0x0002, 0x4000, 0x1000, 0x0400, 0x0100, 0x0040, 0x0010,
         0x0004, 0x0001]
    for x in range((m_array.shape)[0]):
        for y in range((m_array.shape)[1]):
            if x & l[y]:
                m_array[x:x + 1, y:y + 1] = -2
            else:
                m_array[x:x + 1, y:y + 1] = 2

    return m_array


def gen_waveforms(buffer):
    byte_size = len(buffer)
    byte_offset = 0

    multiply_array = build_sign_array()

    # Iterate on building data.
    while byte_offset < byte_size:
        # Read the header...
        (r_len, r_nbr, dri_level, plug_id, r_time, n_subnet, res, dest_plug_id, r_maintype, sr_desc) = \
            unpack_from(datex_hdr_fmt, buffer[byte_offset:])

        # Unpack sub-records info into a list of tuples...
        sub_recs = [unpack_from('<HB', sr_desc[(x * 3):]) for x in range(8)]

        # Inspect for ecg12_wf data.
        if r_maintype == DRI_MT_WAVE and (0, DRI_WF_ECG12) in sub_recs:
            (act_len, status, label, size_in_bytes, status_1, status_2, status_3, status_4, status_5, status_6) = \
                unpack_from(wf_hdr_fmt, buffer[byte_offset + 40:])

            # Ok, we got us a genuine ecg_12 waveform packet
            # Cut a 200mS array for all leads, and an extra one for pace.
            # Pace is a boolean that will be inserted at the correct sample point based on the parse.
            # 100 rows for 200mS (2mS per row).
            leads = np.zeros((900,), dtype=np.int16).reshape((100, 9))
            row = 0

            # Parse samples by looping 20 times..
            for lead in iter_unpack(ecg12_wf_data_fmt, buffer[byte_offset + 60: byte_offset + 1180]):
                a = np.asarray(lead, dtype=np.int16)  # move the tuple to a high spead array
                sign1 = a[40]
                sign2 = a[41]
                pace = a[42] >> 13

                # Get our samples arranged
                b = a[:40].reshape((8, 5)).T
                b[[1, 2], ] = b[[2, 1], ]
                b[[3, 4], ] = b[[4, 3], ]

                # Add correct sign
                b[1, ] *= multiply_array[sign1, 0:8]
                b[2, ] *= multiply_array[sign1, 8:16]
                b[3, ] *= multiply_array[sign2, 0:8]
                b[4, ] *= multiply_array[sign2, 8:16]

                # Finally, add the deltas
                b[1, ] += b[0, ]
                b[2, ] += b[1, ]
                b[3, ] += b[2, ]
                b[4, ] += b[3, ]

                # Insert pace (if present)
                if pace:
                    leads[row + (pace - 1), 8] = 1  # 1 means pace occured

                # Copy to packet
                np.copyto(leads[row:row + 5, 0:8], b)
                row += 5

            yield leads

        # Offset to next sub-record
        byte_offset += r_len


def main(arg_list=None):
    # region Command Line
    ###############################################################################
    # Get command line arguments.
    parser = argparse.ArgumentParser(description="Process the given Datex-Ohmeda Record")
    parser.add_argument('-i', dest='drc_input_file',
                        help='name of Datex-Ohmeda Record to read and process', required=True)
    parser.add_argument('-v', dest='verbose', default=False, action='store_true',
                        help='verbose output flag', required=False)
    parser.add_argument('--version', action='version', help='Print version.',
                        version='%(prog)s Version {version}'.format(version=__version__))

    # Parse the command line arguments
    args = parser.parse_args(arg_list)

    ###############################################################################
    # Test for existence of the LeCroy file.
    if os.path.isfile(args.drc_input_file) is False:
        print('ERROR, ' + args.csv_input_file + ' does not exist')
        print('\n\n')
        parser.print_help()
        return -1

    # endregion

    ###############################################################################
    # Open DRI and read in entire record...
    with open(args.drc_input_file, 'rb') as fh:
        file_data = fh.read()

    if args.verbose is True:
        print('Size of {0:s} is {1:s} bytes'.format(args.drc_input_file, humansize(len(file_data))))

    ###############################################################################
    # Open file and parse.
    output_filename = os.path.splitext(os.path.basename(args.drc_input_file))[0] + '.csv'
    with open(output_filename, 'wb') as fh:
        total_time = 0
        for wf_data in gen_waveforms(file_data):
            np.savetxt(fh, wf_data, fmt='%d', delimiter=',')
            total_time += 1
            if total_time % 100 == 0:
                print('{0:.2f} seconds processed'.format(total_time/100))

    return 0


###############################################################################
if __name__ == '__main__':
    rv = main()
    exit(rv)
