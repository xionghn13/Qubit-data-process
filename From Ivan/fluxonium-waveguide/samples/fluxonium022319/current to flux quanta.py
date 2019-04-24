I0 = 0.9  # mA
hI = 1.52  # mA
I = 2.32
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)