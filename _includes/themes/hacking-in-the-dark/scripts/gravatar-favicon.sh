#!/bin/bash
#
# Original Author: Jason Graham
# Source: http://jason.the-graham.com/2010/10/17/generating_gravatar_favicons_for_jekyll/
#
# This script downloads my gravatar image and uses it to create the
# favicon for the website.

# Replace 'HASH' with the value specific to your email address.
# Also, I pull the biggest version of the image available (256x256),
# but this probably isn't necessary.

# I (akajontsai-devel@yahoo.com) modified this to take a command-line argument

EMAIL=$1
if [ "${EMAIL}" = "" ]
then
    echo "Usage: ./gravatar-favicon.sh EMAIL"
    exit
fi

HASH=`python -c "import hashlib;print hashlib.md5('${EMAIL}').hexdigest()"`

wget http://www.gravatar.com/avatar/${HASH}?s=256 -O /tmp/favicon.jpg

convert /tmp/favicon.jpg ./favicon.png

# below is adapted from
# http://www.brennan.id.au/13-Apache_Web_Server.html#favicon
pngtopnm -mix ./favicon.png > tmp_logo.pnm

rm -f favicon.ico

pnmscale -xsize=32 -ysize=32 ./tmp_logo.pnm > tmp_logo32.ppm
pnmscale -xsize=16 -ysize=16 ./tmp_logo.pnm > tmp_logo16.ppm

pnmquant 256 tmp_logo32.ppm > tmp_logo32x32.ppm
pnmquant 256 tmp_logo16.ppm > tmp_logo16x16.ppm

ppmtowinicon -output ./favicon.ico tmp_logo16x16.ppm tmp_logo32x32.ppm

rm -f tmp_logo* favicon.png
