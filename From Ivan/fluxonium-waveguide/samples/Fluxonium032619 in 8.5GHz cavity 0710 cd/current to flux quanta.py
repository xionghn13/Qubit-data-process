I0 = -82.5e-3  # mA
hI = 2.571  # mA
I = 2.264
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = 8.5
I = flux * hI * 2 + I0
print('I = %.3G' % I)
