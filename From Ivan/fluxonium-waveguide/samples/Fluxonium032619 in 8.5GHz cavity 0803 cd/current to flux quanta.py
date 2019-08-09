I0 = -47.5e-3  # mA
hI = 2.5975  # mA
I = 2.13
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = 8.5
I = flux * hI * 2 + I0
print('I = %.3G' % I)
