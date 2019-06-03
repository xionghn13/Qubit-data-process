I0 = 2.199
hI = 2.469  # mA
I = 5.269
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = -0.0286
I = flux * hI * 2 + I0
print('I = %.3G' % I)
