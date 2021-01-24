var fs = require('fs')
img = fs.readFileSync('test.png')
str = img.toString('base64');