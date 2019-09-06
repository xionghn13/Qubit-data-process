I0 = 32.5e-3  # mA
hI = 2.5975  # mA
I = 1.9
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = 0.425
I = flux * hI * 2 + I0
print('I = %.3G' % I)
