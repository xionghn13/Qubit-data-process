hI = 2.505 - 0.209  # mA
I0 = 4.725 - hI  # mA
I = 4.75
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = .5
I = flux * hI * 2 + I0
print('I = %.5G' % I)
