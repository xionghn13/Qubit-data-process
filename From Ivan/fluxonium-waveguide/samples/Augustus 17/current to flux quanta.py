hI = (1.501 - 0.6675)  # mA
I0 = 0.6675  # mA
I = 1.5031
flux = (I - I0) / (hI * 2) + 0.5 - 0.4964
print('flux = %.5G' % flux)

flux = .1
I = flux * hI * 2 + I0
print('I = %.5G' % I)
