import qrcode
qr = qrcode.make('https://github.com/sielemann/uebertragungsfkt')
qr.save('Bilder/github_qr.png')