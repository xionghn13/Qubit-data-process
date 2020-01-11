#!/usr/bin/env python3
import numpy as np
from copy import copy
from sequence import Sequence
import gates

class CustomSequence(Sequence):
    def generate_sequence(self, config):
        """Sequence for measuring Ramsey as a function of flux using local tunable flux"""

        # get parameters
        dt = config.get('Pulse spacing')
        self.add_gate_to_all(gates.X2p)
        self.add_gate_to_all(gates.Zp, dt=dt)
        self.add_gate_to_all(gates.X2p, dt=dt)


if __name__ == '__main__':
    pass