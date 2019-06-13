I0 = 2.65
hI = 3.7 # mA
I = 4.0845
flux = (I - I0) / (hI * 2)
print('flux = %.3G' % flux)

flux = 0.3
I = flux * hI * 2 + I0
print('I = %.3G' % I)
