I0 = 2.28  # mA
hI = 2.64  # mA
I = 2.129
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = -0.0286
I = flux * hI * 2 + I0
print('I = %.3G' % I)