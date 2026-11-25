# -*- coding: utf-8 -*-
"""Regenerate the images in tdcrpy/assets/ used by the GUI (tdcrpy.gui).

Run from the repository root:  python docs/make_gui_assets.py [bipm_logo.png]

Sources: docs/logo1.png (TDCRPy) and, optionally, BIPM's white header logo
(https://www.bipm.org/o/bipm-theme/images/logo.png), recoloured to the
TDCRPy navy so it reads on the GUI's light header. Without the BIPM file the
existing tdcrpy/assets/bipm_logo.png is left alone.
"""
import os
import sys

from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'tdcrpy', 'assets')
NAVY = (11, 42, 85)


def _transparent_white(im, threshold=246):
    im = im.convert('RGBA')
    px = im.load()
    for y in range(im.height):
        for x in range(im.width):
            r, g, b, a = px[x, y]
            if min(r, g, b) >= threshold:
                px[x, y] = (255, 255, 255, 0)
    return im


def _fit_height(im, h):
    w = round(im.width * h / im.height)
    return im.resize((w, h), Image.LANCZOS)


def tdcrpy_assets():
    logo = _transparent_white(Image.open(os.path.join(ROOT, 'docs',
                                                      'logo1.png')))
    logo = logo.crop(logo.getbbox())
    # 1x and 2x: the window picks the 2x set on a high-DPI screen.
    _fit_height(logo, 56).save(os.path.join(OUT, 'tdcrpy_logo_header.png'))
    _fit_height(logo, 112).save(os.path.join(OUT,
                                             'tdcrpy_logo_header_2x.png'))
    _fit_height(logo, 96).save(os.path.join(OUT, 'tdcrpy_logo_about.png'))
    _fit_height(logo, 192).save(os.path.join(OUT,
                                             'tdcrpy_logo_about_2x.png'))

    # The vial alone, on a square canvas, for the window and taskbar icon.
    left = logo.crop((0, 0, round(logo.width * 0.22), logo.height))
    left = left.crop(left.getbbox())
    side = max(left.size)
    square = Image.new('RGBA', (side, side), (255, 255, 255, 0))
    square.paste(left, ((side - left.width) // 2, (side - left.height) // 2))
    for n in (16, 32, 48, 256):
        square.resize((n, n), Image.LANCZOS).save(
            os.path.join(OUT, f'tdcrpy_icon_{n}.png'))


def bipm_asset(path):
    im = Image.open(path).convert('RGBA')
    px = im.load()
    for y in range(im.height):
        for x in range(im.width):
            a = px[x, y][3]
            px[x, y] = NAVY + (a,)
    im = im.crop(im.getbbox())
    _fit_height(im, 48).save(os.path.join(OUT, 'bipm_logo.png'))
    _fit_height(im, 96).save(os.path.join(OUT, 'bipm_logo_2x.png'))


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    tdcrpy_assets()
    if len(sys.argv) > 1:
        bipm_asset(sys.argv[1])
    print('written to', OUT)
