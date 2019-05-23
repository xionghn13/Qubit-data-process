I0 = 2.351
hI = 3.744  # mA

I = 4.60188
flux = (I - I0) / (hI * 2)
print('flux = %.5G' % flux)

flux = -0.0286
I = flux * hI * 2 + I0
print('I = %.3G' % I)
