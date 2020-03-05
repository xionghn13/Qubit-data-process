import numpy as np
import scipy.interpolate as itp
import matplotlib.pyplot as plt
import SubtractBackgroundFunc as sbf
import QubitSpectrumFunc as qsf
from scipy.optimize import curve_fit
from FunctionLib import T1_curve, rabi_curve, FitTransientTime
import ExtractDataFunc as edf
import h5py

DataPath = 'E:/Projects\Fluxonium\data_process/Fluxonium123118/'
# BackgroundFile = 'one_tone_8.6GHz_to_8.9GHz_-15dBm_4.9mA_10us integration_100Kavg_50KHz step_011319.dat'
Plus50MHzBackgroundFile = 'one_tone_4.09GHz_to_4.3GHz_-15dBm_4.9mA_10us integration_100Kavg_50KHz step_011619.dat'
Plus50MHzBackgroundFile = '011719_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit1.2045GHz_20dBm_4.9_mA_I cos Q sin mod true interleafing_odd readout even ref_avg300k_Rabi1000_duty100000readout3us.h5'
# Plus50MHzBackgroundFile = 'one_tone_8.6GHz_to_8.77GHz_-15dBm_4.7mA_10us integration_100Kavg_50KHz step_012919.dat'
Minus50MHzBackgroundFile = 'one_tone_4.09GHz_to_4.3GHz_-15dBm_4.9mA_10us integration_100Kavg_50KHz step_011419-50MHz.dat'
# Minus50MHzBackgroundFile = 'one_tone_8.6GHz_to_8.77GHz_-15dBm_4.7mA_10us integration_100Kavg_50KHz step_012919.dat'
# RabiFileList = [
#     '011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-5dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi10000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-6dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi10000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-7dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi10000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-8dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi10000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-9dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi10000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-10dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi10000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-11dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-12dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-13dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-14dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-15dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-16dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-17dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-18dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-19dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-29dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-28dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-27dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-26dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-25dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-24dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-23dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-21dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-22dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011219_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-40dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011119_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-30dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi100000_duty200000readout3us.h5',
# 	'011119_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-15dBm_qubit4.104GHz_-30dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# 	'011119_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-30dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg200k_Rabi50000_duty100000readout3us.h5',
# ]
#
RabiFileList = [
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-5dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1200k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-6dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1200k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-7dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1200k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-8dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1200k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-9dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1200k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-10dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1000k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-11dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1000k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-12dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1000k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-13dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1000k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-14dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg1000k_Rabi10000_duty50000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-15dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg800k_Rabi40000_duty100000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-16dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg800k_Rabi40000_duty100000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-17dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg800k_Rabi40000_duty100000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-18dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg800k_Rabi40000_duty100000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-19dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg800k_Rabi40000_duty100000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-20dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg800k_Rabi40000_duty100000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-21dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg800k_Rabi40000_duty100000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-22dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg800k_Rabi40000_duty100000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-23dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg800k_Rabi40000_duty100000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-24dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-25dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-26dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-27dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-28dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-29dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    '011819_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-30dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-35dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    # '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-34dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    # '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-33dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-32dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
    '011919_rabi_CH2(AWG1Vpp)_no pump_readout_4.154GHz__-20dBm_qubit4.104GHz_-31dBm_5.28_mA_I cos Q sin mod true interleafing_odd readout even ref_avg400k_Rabi100000_duty150000readout3us.h5',
]




# RabiFileList = [
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-25dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg100k_Rabi300000_duty800000readout12us.h5',
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-23dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg100k_Rabi300000_duty800000readout12us.h5',
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-21dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg100k_Rabi200000_duty700000readout12us.h5',
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-19dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi200000_duty700000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-17dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi300000_duty800000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-15dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi300000_duty800000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-13dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi300000_duty800000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-11dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi200000_duty700000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-9dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi200000_duty700000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-7dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi200000_duty700000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-5dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi200000_duty700000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-3dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi200000_duty700000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_-1dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi200000_duty700000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_1dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi100000_duty600000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_3dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi30000_duty500000readout12us.h5',
#     '021419_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_5dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi30000_duty500000readout12us.h5',
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_15dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg100k_Rabi1000_duty500000readout12us.h5',
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_13dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi1000_duty500000readout12us.h5',
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_11dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg100k_Rabi5000_duty500000readout12us.h5',
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_9dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg100k_Rabi5000_duty500000readout12us.h5',
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_7dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi8000_duty500000readout12us.h5',
#     '021519_rabi_CH2(AWG1Vpp)_no pump_readout_4.077GHz__-15dBm_qubit4.027GHz_6dBm_0.89_mA_I cos Q sin mod true interleafing_odd readout even ref_avg50k_Rabi8000_duty500000readout12us.h5',
# ]

IQModFreq = 0.05
FitForGamma = True
Gamma_r = 8.72 * 2 * np.pi
FitCorrectedR = False
LogScale = False
FitTransient = True

PhaseSlope = 326.7041108065019
PhaseRefrenceFreq = 4.105

NumFile = len(RabiFileList)
DrivePowerArray = np.zeros([NumFile, ])
# analyze background file
if Plus50MHzBackgroundFile.startswith('one_tone'):
    [Plus50MHzBackFreq, Plus50MHzBackComplex] = edf.readFSweepDat(DataPath + Plus50MHzBackgroundFile)
    Plus50MHzBackPowerStr = Plus50MHzBackgroundFile.split('_')[5][:-3]
    Plus50MHzBackPower = float(Plus50MHzBackPowerStr)

    [Minus50MHzBackFreq, Minus50MHzBackComplex] = edf.readFSweepDat(DataPath + Minus50MHzBackgroundFile)
    Minus50MHzBackPowerStr = Minus50MHzBackgroundFile.split('_')[5][:-3]
    Minus50MHzBackPower = float(Minus50MHzBackPowerStr)
else:
    BackgroundFileStrList = Plus50MHzBackgroundFile.split('_')
    MeasurementType = BackgroundFileStrList[1]
    if MeasurementType == 'rabi' and BackgroundFileStrList[3] != 'no pump':
        MeasurementType = 'pump rabi'
    BackFreqDict = {'t1': 3, 'rabi': 5, 'pump rabi': 6}
    BackFreqInd = BackFreqDict[MeasurementType]
    BackPowerDict = {'t1': 5, 'rabi': 7, 'pump rabi': 8}
    BackPowerInd = BackPowerDict[MeasurementType]
    BackPower = float(BackgroundFileStrList[BackPowerInd][:-3])
    BackLowerFreq = round(float(BackgroundFileStrList[BackFreqInd][:-3]) - IQModFreq, 4)
    BackHigherFreq = round(BackLowerFreq + IQModFreq * 2, 15)
    BackPower = float(BackgroundFileStrList[BackPowerInd][:-3])
    Plus50MHzBackPower = BackPower
    Minus50MHzBackPower = BackPower
    if MeasurementType in ('rabi', 'transient'):
        [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readRabiH5(DataPath + Plus50MHzBackgroundFile)
    elif MeasurementType == 't1':
        [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readT1H5(DataPath + Plus50MHzBackgroundFile)
    Plus50MHzBackFreq = np.array([BackLowerFreq, BackHigherFreq])
    Plus50MHzBackComplex = np.array([ComplexLowerFreq.mean(), ComplexHigherFreq.mean()])
    Minus50MHzBackFreq = Plus50MHzBackFreq
    Minus50MHzBackComplex = Plus50MHzBackComplex

for i, RabiFile in enumerate(RabiFileList):
    RabiFileStrList = RabiFile.split('_')
    MeasurementType = RabiFileStrList[1]
    if MeasurementType == 'rabi' and RabiFileStrList[3] != 'no pump':
        MeasurementType = 'pump rabi'

    ReadoutFreqDict = {'rabi': 5, 'pump rabi': 6}
    ReadoutFreqInd = ReadoutFreqDict[MeasurementType]
    ReadoutPowerDict = {'rabi': 7, 'pump rabi': 8}
    ReadoutPowerInd = ReadoutPowerDict[MeasurementType]
    ReadoutLowerFreq = round(float(RabiFileStrList[ReadoutFreqInd][:-3]) - IQModFreq, 4)
    ReadoutHigherFreq = round(ReadoutLowerFreq + IQModFreq * 2, 15)
    ReadoutPower = float(RabiFileStrList[ReadoutPowerInd][:-3])
    QubitFreqDict = {'rabi': 8}
    QubitFreqInd = QubitFreqDict[MeasurementType]
    QubitFreq = float(RabiFileStrList[QubitFreqInd][5:-3])
    if QubitFreq in (ReadoutLowerFreq, ReadoutHigherFreq):
        MeasurementType = 'transient'
    DrivePowerDict = {'transient': 9}
    DrivePowerInd = DrivePowerDict[MeasurementType]
    DrivePower = float(RabiFileStrList[DrivePowerInd][:-3])

    if MeasurementType in ('rabi', 'transient'):
        [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readRabiH5(DataPath + RabiFile)
    elif MeasurementType == 't1':
        [Time, ComplexLowerFreq, ComplexHigherFreq] = edf.readT1H5(DataPath + RabiFile)

    RComplexLowerFreq = sbf.FPSweepBackgroundCalibrate(ReadoutLowerFreq, ReadoutPower, ComplexLowerFreq,
                                                       Plus50MHzBackFreq,
                                                       Plus50MHzBackComplex, Plus50MHzBackPower)
    RComplexHigherFreq = sbf.FPSweepBackgroundCalibrate(ReadoutHigherFreq, ReadoutPower, ComplexHigherFreq,
                                                        Minus50MHzBackFreq,
                                                        Minus50MHzBackComplex, Minus50MHzBackPower)
    ComplexLowerFreqNormalized = ComplexLowerFreq * 10 ** (- ReadoutPower / 20)
    ComplexHigherFreqNormalized = ComplexHigherFreq * 10 ** (- ReadoutPower / 20)

    if FitCorrectedR:
        y_data = np.array(np.real(RComplexLowerFreq / RComplexHigherFreq), dtype='float64')
    else:
        y_data = np.array(RComplexLowerFreq.real, dtype='float64')
    x_data = np.array(Time, dtype='float64')
    if MeasurementType in ('t1', 'transient'):
        B_guess = y_data[-1]
        A_guess = y_data[0].real - B_guess
        T1_guess = x_data[-1] / 2
        bounds = (
            (-2, 1, -1),
            (2, 1e6, 1)
        )
        opt, cov = curve_fit(T1_curve, x_data, y_data, p0=[A_guess, T1_guess, B_guess], maxfev=30000)
        A_fit, T1_fit, B_fit = opt
        FitR = T1_curve(Time, A_fit, T1_fit, B_fit)
        ParamList = ['A', 'Decay time/ns', 'B']
    elif MeasurementType in ('rabi'):
        B_guess = y_data.mean()
        A_guess = y_data[0] - B_guess
        T1_guess = x_data[-1]
        MaxInd = y_data.argmax()
        MinInd = y_data.argmin()
        Tpi_guess = np.abs(x_data[MaxInd] - x_data[MinInd])
        phi0_guess = 0
        guess = ([A_guess, T1_guess, B_guess, Tpi_guess, phi0_guess])
        bounds = (
            (-2, 1, -1, 1, - np.pi / 2),
            (2, np.inf, 1, np.inf, np.pi / 2)
        )
        opt, cov = curve_fit(rabi_curve, x_data, y_data, p0=guess, bounds=bounds)
        A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit = opt
        FitR = rabi_curve(Time, A_fit, T1_fit, B_fit, Tpi_fit, phi0_fit)
        ParamList = ['A', 'Decay time/ns', 'B', 'Tpi', 'phi0']
    if i == 0:
        TimeList = []
        RComplexLowerFreqList = []
        RComplexHigherFreqList = []
        FitRList = []
        OptMatrix = np.zeros([len(opt), NumFile])
        ErrMatrix = np.zeros([len(cov), NumFile])
    TimeList.append(Time)
    RComplexLowerFreqList.append(RComplexLowerFreq)
    RComplexHigherFreqList.append(RComplexHigherFreq)
    FitRList.append(FitR)
    ErrMatrix[:, i] = np.sqrt(cov.diagonal())
    for i_p, par in enumerate(opt):
        if ErrMatrix[i_p, i] > par or opt[1] > 5 * Time[-1]:
            OptMatrix[i_p, i] = np.nan
        else:
            OptMatrix[i_p, i] = par

    DrivePowerArray[i] = DrivePower


#sort array
SortInd = np.argsort(DrivePowerArray)
DrivePowerArray = np.sort(DrivePowerArray)
TimeList = [TimeList[i] for i in SortInd]
RComplexLowerFreqList = [RComplexLowerFreqList[i] for i in SortInd]
RComplexHigherFreqList = [RComplexHigherFreqList[i] for i in SortInd]
FitRList = [FitRList[i] for i in SortInd]
OptMatrix = OptMatrix[:, SortInd]
ErrMatrix = ErrMatrix[:, SortInd]

if FitTransient:
    opt, cov, FitTime = FitTransientTime(DrivePowerArray, Gamma_r, OptMatrix)
    Gamma_in_fit, Gamma_out_fit, pow_ratio_fit = opt
    TransientOpt = opt
    TransientErr = np.sqrt(cov.diagonal())
    TransientParamList = ['Gamma_in', 'Gamma_out', 'pow_ratio']

limit = 1.7

fig, ax = plt.subplots()
for i in range(NumFile):
    plt.plot(np.real(RComplexLowerFreqList[i]), np.imag(RComplexLowerFreqList[i]))
    plt.plot(np.real(RComplexHigherFreqList[i]), np.imag(RComplexHigherFreqList[i]))
plt.plot([-2, 2], [0, 0], '--')
plt.plot([1], [0], 'ro')
plt.xlabel('Re', fontsize='x-large')
plt.ylabel('Im', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.xlim(-limit, limit)
plt.ylim(-limit, limit)
ax.set_aspect('equal')

fig, ax = plt.subplots()
for i in range(NumFile):
    plt.plot(np.real(RComplexLowerFreqList[i] / RComplexHigherFreqList[i]),
             np.imag(RComplexLowerFreqList[i] / RComplexHigherFreqList[i]))
plt.plot([-2, 2], [0, 0], '--')
plt.plot([1], [0], 'ro')
plt.xlabel('Re', fontsize='x-large')
plt.ylabel('Im', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
plt.xlim(-limit, limit)
plt.ylim(-limit, limit)
ax.set_aspect('equal')

fig, ax = plt.subplots()
for i in range(NumFile):
    plt.plot(TimeList[i], np.real(RComplexLowerFreqList[i]), 'o')
    if FitCorrectedR:
        plt.plot(TimeList[i], np.real(RComplexHigherFreqList[i]), 'o')
    if not FitCorrectedR:
        plt.plot(TimeList[i], FitRList[i])
plt.xlabel('Time/ns', fontsize='x-large')
plt.ylabel('Re', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
for i in range(NumFile):
    plt.plot(TimeList[i], np.real(RComplexLowerFreqList[i] / RComplexHigherFreqList[i]), 'o')
    if FitCorrectedR:
        plt.plot(TimeList[i], FitRList[i])
plt.xlabel('Time/ns', fontsize='x-large')
plt.ylabel('Re', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()

fig, ax = plt.subplots()
plotInd = 1
# plt.plot(DrivePowerArray, OptMatrix[plotInd, :]/1000, 'o')
ax.errorbar(DrivePowerArray, OptMatrix[plotInd, :] / 1000, yerr=ErrMatrix[plotInd, :] / 1000, fmt='o')
if FitTransient:
    plt.plot(DrivePowerArray, FitTime)
plt.xlabel('Power/dBm', fontsize='x-large')
plt.ylabel('Decay time/us', fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
if FitTransient:
    plt.title('Gamma_in=%.3G$\pm$%.3GMHz, Gamma_out=%.3G$\pm$%.3GMHz,\n pow_ratio=%.3G$\pm$%.3G' % (
        Gamma_in_fit, TransientErr[0], Gamma_out_fit, TransientErr[1], pow_ratio_fit, TransientErr[2]))
plt.tight_layout()
if LogScale:
    ax.set_yscale('log')
# plt.ylim(0, 40000)


fig, ax = plt.subplots()
plotInd = 2
# plt.plot(DrivePowerArray, OptMatrix[plotInd, :]/1000, 'o')
ax.errorbar(DrivePowerArray, OptMatrix[plotInd, :], yerr=ErrMatrix[plotInd, :], fmt='o')
plt.xlabel('Power/dBm', fontsize='x-large')
plt.ylabel(ParamList[plotInd], fontsize='x-large')
plt.tick_params(axis='both', which='major', labelsize='x-large')
plt.tight_layout()
# plt.ylim(0, 20000)

plt.show()
