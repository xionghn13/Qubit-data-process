I0 = 2.351
hI = 3.744  # mA

I = 4.59754
flux = (I - I0) / (hI * 2)
print('flux = %.5G' % flux)

flux = 2.5
I = flux * hI * 2 + I0
print('I = %.3G' % I)
