"""Module with color and colormap definitions.

Colors are from colorbrewer2.org, except for 'white' and 'black'. They were
obtained with these settings:
Number of data classes = 7
Nature of data = qualitative
Color scheme = Set1

Colormaps were exported from Paraview 4.3.1

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


def white():
    """The color white.

    HEX = #FFFFFF
    RGB = 255, 255, 255
    CMYK = 0, 0, 0, 0

    """
    rgb = [255, 255, 255]
    return [i / 255. for i in rgb]



def black():
    """The color white.

    HEX = #000000
    RGB = 0, 0, 0
    CMYK = 0, 0, 0, 1

    """
    rgb = [0, 0, 0]
    return [i / 255. for i in rgb]


def BuRd():
    """Blue to red colormap.

    BuRd in CIELAB color space.

    """
    colors = [[-1, 0.0196078, 0.188235, 0.380392],
              [-0.87451, 0.0885023, 0.321111, 0.564935],
              [-0.74902, 0.163394, 0.444984, 0.697505],
              [-0.623529, 0.247059, 0.555703, 0.754101],
              [-0.498039, 0.420691, 0.676432, 0.818692],
              [-0.372549, 0.606455, 0.789776, 0.880278],
              [-0.247059, 0.761471, 0.868513, 0.924559],
              [-0.121569, 0.87805, 0.925719, 0.951949],
              [0.00392157, 0.969085, 0.966476, 0.964935],
              [0.129412, 0.983856, 0.897581, 0.84683],
              [0.254902, 0.982467, 0.800687, 0.706111],
              [0.380392, 0.960327, 0.667826, 0.536339],
              [0.505882, 0.894575, 0.503807, 0.399771],
              [0.631373, 0.817075, 0.332174, 0.281041],
              [0.756863, 0.728496, 0.155016, 0.197391],
              [0.882353, 0.576928, 0.0553597, 0.149248],
              [1, 0.403922, 0, 0.121569]]
    return colors


def coolwarm():
    """coolwarm colormap.

    coolwarm in CIELAB color space

    """
    colors = [[-1, 0.229801, 0.298711, 0.753689]
              [-0.875, 0.303868, 0.406531, 0.844953]
              [-0.75, 0.383017, 0.509422, 0.917388]
              [-0.625, 0.466667, 0.604562, 0.968154]
              [-0.5, 0.552956, 0.68893, 0.995377]
              [-0.375, 0.63917, 0.759594, 0.998154]
              [-0.25, 0.722194, 0.813947, 0.976577]
              [-0.125, 0.798688, 0.84979, 0.931685]
              [0, 0.8654, 0.865415, 0.8654]
              [0.125, 0.924132, 0.82739, 0.774502]
              [0.25, 0.958846, 0.769772, 0.678004]
              [0.375, 0.969955, 0.69427, 0.57937]
              [0.5, 0.958007, 0.602838, 0.481773]
              [0.625, 0.923949, 0.497307, 0.387976]
              [0.75, 0.869184, 0.378317, 0.300267]
              [0.875, 0.795636, 0.241291, 0.220523]
              [1, 0.705669, 0.0155489, 0.15024]]
