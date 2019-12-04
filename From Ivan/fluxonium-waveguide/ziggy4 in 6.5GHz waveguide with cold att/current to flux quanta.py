hI = 2.505 - 0.209  # mA
I0 = 0.13 - hI  # mA
I = I0 + 0.031
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = 1.5
I = flux * hI * 2 + I0
print('I = %.5G' % I)
