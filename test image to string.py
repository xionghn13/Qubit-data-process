import pytesseract

from PIL import Image

value = Image.open("test.png")

text = pytesseract.image_to_string(value, config='--tessdata-dir "C://Program Files//Tesseract-OCR//tessdata"')

print("text present in images:", text)