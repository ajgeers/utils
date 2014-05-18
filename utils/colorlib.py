"""Module with color definitions from colorbrewer2.org

The colors were obtained with these settings:
Number of data classes = 7
Nature of data = qualitative
Color scheme = Set1

"""

def red():
    """The color red.

    HEX = #e41a1c
    RGB = 228, 26, 28
    CMYK = 10, 90, 80, 0

    """
    rgb = [228, 26, 28]
    return [i / 255. for i in rgb]


def blue():
    """The color blue.

    HEX = #377eb8
    RGB = 55, 126, 184
    CMYK = 80, 30, 0, 0

    """
    rgb = [55, 126, 184]
    return [i / 255. for i in rgb]


def green():
    """The color green.

    HEX = #4daf4a
    RGB = 77, 175, 74
    CMYK = 70, 0, 80, 0

    """
    rgb = [77, 175, 74]
    return [i / 255. for i in rgb]


def purple():
    """The color purple.

    HEX = #984ea3
    RGB = 152, 78, 163
    CMYK = 40, 65, 0, 0

    """
    rgb = [152, 78, 163]
    return [i / 255. for i in rgb]


def orange():
    """The color orange.

    HEX = #ff7f00
    RGB = 255, 127, 0
    CMYK = 0, 50, 100, 0

    """
    rgb = [255, 127, 0]
    return [i / 255. for i in rgb]


def yellow():
    """The color yellow.

    HEX = #ffff33
    RGB = 255, 255, 51
    CMYK = 0, 0, 80, 0

    """
    rgb = [255, 255, 51]
    return [i / 255. for i in rgb]


def brown():
    """The color brown.

    HEX = #a65628
    RGB = 166, 86, 40
    CMYK = 35, 60, 80, 0

    """
    rgb = [166, 86, 40]
    return [i / 255. for i in rgb]
