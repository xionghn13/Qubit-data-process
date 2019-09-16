I0 = 3.94e-3  # mA
hI = 2.61806  # mA
I = 1.
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = 0.425
I = flux * hI * 2 + I0
print('I = %.3G' % I)
