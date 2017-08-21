###############################################################################
# Python script used to read in Datex-Ohmeda Record file and parse out
# ECG12 sub-records.

###############################################################################

# region Import region
import os
import argparse  # command line parser
from struct import *
from datetime import datetime
import tags
# endregion

# region Globals region

__version__ = '1.0'  # version of script

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

# endregion


def parse_file(file_handle, main_record_type, sub_record_type):
    offset = 0     # we need to keep track of where we are in the file
    while True:
        header_data = file_handle.read(DRC_HEADER_SIZE)
        if len(header_data) == 0:
            break
        (r_len, r_nbr, dri_level, plug_id, r_time, n_subnet, res, dest_plug_id, r_maintype, sr_desc) = \
            unpack_from('<HBBHLBBHH24s', header_data)
        # Construct Datex - Record header
        drh = \
        {
            'r_len': r_len,
            'r_nbr': r_nbr,
            'dri_level': dri_level,
            'plug_id': plug_id,
            'r_time': datetime.fromtimestamp(r_time).strftime('%Y-%m-%d %H:%M:%S'),
            'n_subnet': n_subnet,
            'res': res,
            'dest_plug_id': dest_plug_id,
            'r_maintype': r_maintype,
            'sr_desc': sr_desc
        }

        # Unpack sub-records info into a list of tuples...
        sub_recs = [unpack_from('<HB', sr_desc[(x*3):]) for x in range(8)]

        # ECG waveform?
        if r_maintype == main_record_type:
            # Hunt down the sub-record the caller wants...
            if (0, sub_record_type) in sub_recs:
                # We got us a ECG sub-record to parse!
                file_handle.seek(offset+60)     # offset directly to waveforms
                wf_data = file_handle.read(1120)    # read entire waveform, 5 samples for each of the 8 leads

                for x in range(0, 1120, 56):
                    (Sign1, Sign2, Pacer, QRS) = unpack_from('<40xHHHH', wf_data[x:])

                    # Calc lead I...
                    (s1, s2, s3) = unpack_from('<hHH', wf_data[x:])
                    di2 = (s2 >> 8)
                    if Sign1 & 0x8000:
                        di2 = -di2
                    di3 = (s2 & 0xff)
                    if Sign1 & 0x4000:
                        di3 = -di3
                    di4 = (s3 >> 8)
                    if Sign2 & 0x8000:
                        di4 = -di4
                    di5 = (s3 & 0xff)
                    if Sign2 & 0x4000:
                        di5 = -di5
                    lead_i = [None] * 5
                    lead_i[0] = s1
                    lead_i[1] = lead_i[0] + di2
                    lead_i[2] = lead_i[1] + di3
                    lead_i[3] = lead_i[2] + di4
                    lead_i[4] = lead_i[3] + di5
                    yield lead_i   # send back

        # Advance to next sub-record...
        offset += r_len
        file_handle.seek(offset)


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
    # Open file and parse.
    with open(args.drc_input_file, 'rb') as f:
        for lead_i in parse_file(f, DRI_MT_WAVE, DRI_WF_ECG12):
            for s in lead_i:
                print(s)


###############################################################################
if __name__ == '__main__':
    rv = main()
    exit(rv)
