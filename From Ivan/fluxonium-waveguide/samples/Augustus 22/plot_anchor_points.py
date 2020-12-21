import os
import sys
import numpy as np
import warnings
from matplotlib import pyplot as plt

import Labber

if __file__ in [f for f in os.listdir('.') if os.path.isfile(f)]:
    SCRIPT_PATH = os.path.dirname(os.getcwd())
else:
    SCRIPT_PATH = os.path.dirname(__file__)
LOCAL_PATH = os.path.join(SCRIPT_PATH.rsplit('fluxonium-waveguide', 1)[0],
                                             'fluxonium-waveguide')
if LOCAL_PATH not in sys.path:
    sys.path.append(LOCAL_PATH)

from local_settings import samples_path
import plot_spectrum
import anchor_points


def plot_anchor_points(points, offset=.1):
    if isinstance(points, str):
        points = anchor_points.data[points]
    phi_ext = []
    freqs = []
    labels = []
    plot_labels = True
    for point in points:
        phi_ext.append(point['external flux quanta'])
        freqs.append(point['frequency'])
        if 'label' in point:
            labels.append(point['label'])
        else:
            labels.append('%d-%d' %point['transition'])
    
    plt.scatter(phi_ext, freqs, color='k')
    if plot_labels:
        for idx in range(len(phi_ext)):
            plt.gca().text(phi_ext[idx], freqs[idx] + offset, labels[idx],
                   color='k', fontsize=12,
                   horizontalalignment='center',
                   bbox=dict(boxstyle='round', fc='w', alpha=.5))


def show_labels(points, offset=.1):
    if isinstance(points, str):
        points = anchor_points.data[points]
    phi_ext = []
    freqs = []
    labels = []
    for point in points:
        phi_ext.append(point['external flux quanta'])
        freqs.append(point['frequency'])
        if 'label' in point:
            labels.append(point['label'])
        else:
            labels.append('%d-%d' %point['transition'])

    for idx in range(len(phi_ext)):
        plt.gca().text(phi_ext[idx], freqs[idx] + offset, labels[idx],
               color='k', fontsize=14,
               horizontalalignment='center',
               bbox=dict(boxstyle='round', fc='w', alpha=.5))


if __name__ == '__main__':
    plot_spectrum.plot_spectrum()
    plot_anchor_points('data1')
    plot_spectrum.label_axes(ylim=[0, 10])
    plt.show()
