import ReaderSlice
import ReaderRay
import matplotlib.pyplot as plt

import AnalyticalDualMax

#The location of output files to be read
def PATH(test):
    return f'../test/{test}/output'

#The test directories to read from
TESTS = [
    'SR_RES_32_4',
    'SR_RES_64_8',
    'SR_RES_128_16',
    'SR_RES_256_32'
]

#The dump function name to be read.
#For example 'a_p_slice_512_32_00000000' will be read by 'a_p'
READ_VALUE = 'a_p'

#The type of dump file to read
#Either 'slice' or 'rays'
FILE_TYPE = 'rays'

#Limits of the plot
XLIMITS = [5, 11]
YLIMITS = [-0.5, 0.5]

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
    AnalyticalDualMax.multipole(1, 0, 0)
]

#A list of function(s) that takes arguments (numerical, analytical)
#Only works if there is an analytical function
#Not required if PLOT_AUXILIARY is set to False
AUXILIARY = [
    lambda num, an: num - an,
    lambda num, an: 16 * (num - an),
    lambda num, an: 16 * 16 * (num - an),
    lambda num, an: 16 * 16 * 16 * (num - an)
]

#The horizontal axis of the plot to be created
#HORIZONTAL: 'r', 'theta', or 't'
#HORIZONTAL_GRID: range of gridpoint(s) or None to get all available
HORIZONTAL = 'r'
HORIZONTAL_GRID = None

#Plot multiple different lines in the same plot
#SERIES: 'r', 'theta', or 't'
#SERIES_GRID list containing one or more grid step(s)
SERIES = 'theta'
SERIES_GRID = [1]

#The component that can be animated
#MOVIE: 'r', 'theta', or 't'
#MOVIE_GRID: list containing animation gridpoint(s) or None to get all available
MOVIE = 't'
MOVIE_GRID = [1]

#Scaling the time step for simulations that run at different resolutions
MOVIE_SCALAR = [
    1,
    4,
    16,
    64
]

#initialize list
movie_length = [0 for _ in range(len(TESTS))]

if FILE_TYPE.lower() == 'slice':
    readers = [
        ReaderSlice.SliceReader(PATH(TESTS[i]), READ_VALUE, HORIZONTAL, HORIZONTAL_GRID, SERIES, SERIES_GRID, 
                                MOVIE, [val * MOVIE_SCALAR[i] for val in MOVIE_GRID]) for i in range(len(TESTS))
    ]
elif FILE_TYPE.lower() == 'rays':
     readers = [
        ReaderRay.RayReader(PATH(TESTS[i]), READ_VALUE, HORIZONTAL, HORIZONTAL_GRID, SERIES, SERIES_GRID, 
                                MOVIE, [val * MOVIE_SCALAR[i] for val in MOVIE_GRID]) for i in range(len(TESTS))
    ] 
else:
    print(f"File type {FILE_TYPE} not recognized as 'slice'or 'rays'.")
    raise TypeError


for i, reader in enumerate(readers):
    reader.set_analytical(ANALYTICAL[i])
    reader.set_auxiliary(AUXILIARY[i])

    reader.read()
    reader.set_xlimits(XLIMITS)
    reader.set_ylimits(YLIMITS)
    reader.set_legend_prefix(TESTS[i])

    movie_length[i] = len(reader.z_grid)


if all(length == movie_length[0] for length in movie_length):
    print("All movie steps length match! Plotting . . .")
    _, ax = plt.subplots(1)

    for movie_iteration in range(movie_length[0]):
        ax.clear()
        for i in range(len(TESTS)):
            readers[i].plot_z(movie_iteration, ax, PLOT_NUMERICAL, (PLOT_ANALYTICAL and ANALYTICAL != None), (PLOT_AUXILIARY and AUXILIARY != None))
        plt.pause(0.01)
    plt.show()
else:
    print("Movie length mismatch !!!!!")

del readers