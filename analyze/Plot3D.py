import ReaderSlice
import ReaderRay
import math
import matplotlib.pyplot as plt

import AnalyticalDualMax

#The location of output files to be read
def PATH(test):
    return f'../test/{test}/output'

#The test directories to read from
TESTS = [
    # 'Conv_OD_256_32',
    # 'Conv_OD_384_48',
    # 'Conv_OD_512_64',
    # 'Conv_OD_768_96',
    # 'Conv_OD_1024_128',
    'New_Grid'
]

#The dump function name to be read.
#For example 'a_p_slice_512_32_00000000' will be read by 'a_p'
READ_VALUE = 'a_p'

#The type of dump file to read
#Either 'slice' or 'rays'
FILE_TYPE = 'rays'

#Limits of the plot
XLIMITS = [30,50]
YLIMITS = None

#Whether or not to plut each function
PLOT_NUMERICAL = False
PLOT_ANALYTICAL = False
PLOT_AUXILIARY = True

#A list of function(s) for each test the takes arguments (r, theta, time)
#Not required if PLOT_ANALYTICAL is set to False
ANALYTICAL = [
    AnalyticalDualMax.multipole(1, 0, 0),
    AnalyticalDualMax.multipole(1, 0, 0),
    AnalyticalDualMax.multipole(1, 0, 0),
    AnalyticalDualMax.multipole(1, 0, 0),
    AnalyticalDualMax.multipole(1, 0, 0)
]

#A list of function(s) that takes arguments (numerical, analytical)
#Only works if there is an analytical function
#Not required if PLOT_AUXILIARY is set to False
AUXILIARY = [
    lambda num, an: (num - an),
    lambda num, an: 5.0625 * (num - an),
    lambda num, an: 16 * (num - an),
    lambda num, an: 5.0625 * 16 * (num - an),
    lambda num, an: 16 * 16 * (num - an),

]

# AUXILIARY = [
#     lambda num, an: num - an,
#     lambda num, an: (num - an),
#     lambda num, an: (num - an),
#     lambda num, an: (num - an),
#     lambda num, an: (num - an)
# ]

#The horizontal axis of the plot to be created
#HORIZONTAL: 'r', 'theta', or 't'
#HORIZONTAL_GRID: range of gridpoint(s) or None to get all available
HORIZONTAL = 'r'
HORIZONTAL_GRID = range(0,256)

#Plot multiple different lines in the same plot
#SERIES: 'r', 'theta', or 't'
#SERIES_GRID list containing one or more grid step(s)
SERIES = 'theta'
SERIES_GRID = [1]

#The component that can be animated
#MOVIE: 'r', 'theta', or 't'
#MOVIE_GRID: list containing animation gridpoint(s) or None to get all available
MOVIE = 't'
MOVIE_GRID = range(500,700)
MOVIE_PLOT = None

#initialize list
movie_length = [0 for _ in range(len(TESTS))]

if FILE_TYPE.lower() == 'slice':
    readers = [
        ReaderSlice.SliceReader(PATH(TESTS[i]), READ_VALUE, HORIZONTAL, HORIZONTAL_GRID, SERIES, SERIES_GRID, 
                                MOVIE, MOVIE_GRID) for i in range(len(TESTS))
    ]

elif FILE_TYPE.lower() == 'rays':
     readers = [
        ReaderRay.RayReader(PATH(TESTS[i]), READ_VALUE, HORIZONTAL, HORIZONTAL_GRID, SERIES, SERIES_GRID, 
                                MOVIE, MOVIE_GRID) for i in range(len(TESTS))
    ] 
else:
    print(f"File type {FILE_TYPE} not recognized as 'slice'or 'rays'.")
    raise TypeError

min_resolution = [min([readers[i].resolution[j] for i in range(len(readers))]) for j in range(3)]

for i, reader in enumerate(readers):
    reader.set_analytical(ANALYTICAL[i])
    reader.set_auxiliary(AUXILIARY[i])
    reader.set_xlimits(XLIMITS)
    reader.set_ylimits(YLIMITS)
    reader.set_legend_prefix(TESTS[i])

    reader.rescale_x(min_resolution)
    reader.fill_x()

    reader.read()

    # reader.abs_integrate_plot(reader.aux_plot, 'x')


    if MOVIE_PLOT == None:
        movie_length[i] = len(reader.z_grid)
    else:
        movie_length[i] = len(MOVIE_PLOT)


# maximum = max([readers[i].max() for i in range(len(readers))]) #type: ignore
# for reader in readers: 
#     factor = maximum / reader.max()
#     reader.scale(factor)
#     print(f'Factor: {factor}')


if all(length == movie_length[0] for length in movie_length):
    print("All movie steps length match! Plotting . . .")

    if MOVIE_PLOT == None:
        MOVIE_PLOT = range(movie_length[0])

    fig, ax = plt.subplots(1)

    plot_active = True
    def on_close(event):
        global plot_active
        print("Figure Closed")
        plot_active = False
    fig.canvas.mpl_connect('close_event', on_close)

    for movie_iteration in MOVIE_PLOT:
        if plot_active:
            ax.clear()
            for i in range(len(TESTS)):
                readers[i].plot_z(movie_iteration, ax, PLOT_NUMERICAL, (PLOT_ANALYTICAL and ANALYTICAL != None), (PLOT_AUXILIARY and AUXILIARY != None))
            ax.set_xlim(XLIMITS)
            ax.set_ylim(YLIMITS)
            plt.legend()
            plt.pause(0.01)
        else:
            break
    if plot_active:
        plt.show()
else:
    print("Movie length mismatch !!!!!")

print("Deleting Readers")
del readers