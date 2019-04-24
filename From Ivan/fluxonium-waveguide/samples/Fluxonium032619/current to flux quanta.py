I0 = 2.28  # mA
hI = 2.64  # mA
I = 3.515
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = 0.52
I = flux * hI * 2 + 2.28
print('I = %.3G' % I)