data = {
    'data1': [
        {'transition': (0, 2),  # 00-01
         'external flux quanta': 0.9585,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 5.959  # GHz
         },
        # {'transition': (0, 1),  # 00-10
        #  'external flux quanta': 0.364,  # \Phi_\mathrm{ext}/\Phi_0
        #  'frequency': 3.09  # GHz
        #  },
        # {'transition': (0, 1),  # 00-10
        #  'external flux quanta': 0.438,  # \Phi_\mathrm{ext}/\Phi_0
        #  'frequency': 1.59  # GHz
        #  },
        # {'transition': (0, 1),  # 00-10
        #  'external flux quanta': 0.524,  # \Phi_\mathrm{ext}/\Phi_0
        #  'frequency': 0.892  # GHz
        #  },
        # {'transition': (0, 1),  # 00-10
        #  'external flux quanta': 0.5923,  # \Phi_\mathrm{ext}/\Phi_0
        #  'frequency': 2.22  # GHz
        #  },
        {'transition': (0, 1),  # 00-10
         'external flux quanta': 0.8272,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 4.969  # GHz
         },
        {'transition': (0, 1),  # 00-10
         'external flux quanta': 0.9954,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 5.212  # GHz
         },
        # half flux quanta
        {'transition': (0, 1),  # 00-10
         'external flux quanta': 0.5,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 0.7065  # GHz
         },
        {'transition': (0, 2),  # 00-01
         'external flux quanta': 0.5,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 1.310  # GHz
         },
        {'transition': (0, 3),  # 00-11
         'external flux quanta': 0.5,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 1.310 + 0.7065 - 2.1e-3  # GHz
         },
        {'transition': (3, 6),  # 11-12
         'external flux quanta': 0.5,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 3.205  # GHz
         },
        {'transition': (2, 5),  # 01-02
         'external flux quanta': 0.5,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 3.275  # GHz
         },
        {'transition': (1, 4),  # 10-20
         'external flux quanta': 0.5,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 3.378  # GHz
         },
        {'transition': (3, 7),  # 11-21
         'external flux quanta': 0.5,  # \Phi_\mathrm{ext}/\Phi_0
         'frequency': 3.448  # GHz
         },

        # 0:00, 1:10, 2:01, 3:11, 4:20, 5:02, 6:12, 7:21
        # 11 - 12 frequency: 3.205 GHz, 00 - 12: 5.2194
        # 01 - 02 frequency: 3.275 GHz, 00 - 02: 4.585
        # 10 - 20 frequency: 3.378 GHz, 00 - 20: 4.0845
        # 11 - 21 frequency: 3.448 GHz, 00 - 21: 5.4624

# 11->21: 3.448 GHz, 81 MHz detuning Rabi rate at CZ power would be 190 MHz
# 10->20: 3.378 GHz, 150 MHz detuning Rabi rate at CZ power would be 134 MHz
# 01->02: 3.275 GHz, 253 MHz detuning Rabi rate at CZ power would be 144 MHz
# 11->12: 3.205 GHz, 323 MHz detuning Rabi rate at CZ power would be 75 MHz
# 00-10: 706.5 MHz, 00-01: 1310 MHz
# 11->21: T1 500 ns, T2ramsey 1.2 us, T2echo=1.2 us
# 10->20: T1 400 ns, T2ramsey 475 ns, T2echo=0.8 us
# 01->02: T1 1.2 us, T2ramsey 1.2 us, T2echo=1.7 us
# 11->12: T1 700 ns, T2ramsey 1.2 us, T2echo=1.3 us
#         We have a sigma=12 gaussian part with chop=3, so the gaussian width is 36 ns. We optimized drag for the pulse (slightly more than 3 for the drag coefficient).

],

}
