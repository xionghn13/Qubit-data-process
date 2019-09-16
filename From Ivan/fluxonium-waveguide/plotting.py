import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.collections import LineCollection


colorsC = {0: to_rgba('C7'),
           1: to_rgba('C0'),
           2: to_rgba('C3'),
           3: to_rgba('C4'),
           4: to_rgba('C1'),
           5: to_rgba('C2'),
           6: to_rgba('C6'),
           7: to_rgba('C5'),
           8: to_rgba('C8'),
           9: to_rgba('C9'),
          10: to_rgba('r'),
          11: to_rgba('g'),
          12: to_rgba('m'),
          13: to_rgba('b'),
          14: to_rgba('c'),
          15: to_rgba('k')}
          
colors = { 0: to_rgba('C7'),
           1: to_rgba('b'),
           2: to_rgba('r'),
           3: to_rgba('g'),
           4: to_rgba('m'),
           5: to_rgba('c'),
           6: to_rgba('C0'),
           7: to_rgba('C3'),
           8: to_rgba('C4'),
           9: to_rgba('C1'),
          10: to_rgba('C2'),
          11: to_rgba('C6'),
          12: to_rgba('C5'),
          13: to_rgba('C8'),
          14: to_rgba('C9')}
            
def compute_colors(weights):
    clrs = []
    for k in range(weights.shape[0]):
        clr = np.zeros(4)
        for idx in range(weights.shape[1]):
            clr += weights[k,idx] * np.array(colors[idx % len(colors)])
        clrs.append(tuple(np.clip(clr, 0., 1.)))
    return clrs
    
def colored_lines(x, levels, weights, initial, final, num_photons=1):
    lines = []
    if initial == 0 and num_photons == 1:
        width = 1.5
        style = '-'
    elif initial == 1 and num_photons == 1:
        width = 1.5
        style = '--'
    else:
        width = 1.0
        style = ':'
    for idx in range(initial+1, final):
        freqs = (levels[:,idx] - levels[:,initial]) / float(num_photons)
        points = np.array([x, freqs]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lines.append(LineCollection(segments,
            colors=compute_colors(weights[:,idx]),
            linewidths=width, linestyles=style))
    return lines

def correct_grid(xx, yy):
    xx = xx.T
    yy = yy.T
    dxx = np.diff(xx, axis=1)
    # print(dxx)
    xx = .5 * (np.concatenate([np.array([xx[:,0]-dxx[:,0]]).T, xx],
                               axis=1)
             + np.concatenate([xx, np.array([xx[:,-1]+dxx[:,-1]]).T],
                               axis=1))
    xx = np.concatenate([np.array([xx[0]]), xx], axis=0)

    dyy = np.diff(yy, axis=0)
    yy = .5 * (np.concatenate([np.array([yy[0]-dyy[0]]), yy], axis=0)
             + np.concatenate([yy, np.array([yy[-1]+dyy[-1]])], axis=0))
    yy = np.concatenate([yy, np.array([yy[:,-1]]).T], axis=1)

    xx = xx.T
    yy = yy.T
    return xx, yy
