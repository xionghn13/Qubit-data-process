I0 = - 0.102  # mA
hI = 2.528 + 0.102  # mA
I = 1.
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = 4.5
I = flux * hI * 2 + I0
print('I = %.5G' % I)
