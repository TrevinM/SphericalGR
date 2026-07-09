import AnalyticalDualMax
import ReaderRay
OUTPUT_PATH = '../test/SR_RES_256_32/output'


#The dump function name to be read as a slice file.
#For example 'a_p_slice_512_32_00000000' will be read by 'a_p'
READ_VALUE = 'a_p'

#Vertical axis limits of the plot
LIMITS = [-1, 1]

#A function the takes arguments (r, theta, time)
#Set it as None to not compare to an analytical function
ANALYTICAL = AnalyticalDualMax.multipole(1,0,0,0)

#A function that takes arguments (numerical, analytical)
#Set it as None to not use an auxiliary function
#Only works if there is an analytical function
AUXILIARY = lambda num, an: num - an

#The horizontal axis of the plot to be created
#HORIZONTAL: 'r', 'theta', or 't'
#HORIZONTAL_GRID: range of gridpoints or None to get all available
HORIZONTAL = 'r'
HORIZONTAL_GRID = None

#Plot multiple different lines in the same plot
#SERIES: 'r', 'theta', or 't'
#SERIES_GRID list containing one or more grid step(s)
SERIES = 'theta'
SERIES_GRID = [1]

#The component that can be animated
#MOVIE: 'r', 'theta', or 't'
#MOVIE_GRID: list containing animation values or None to get all available
MOVIE = 't'
MOVIE_GRID = None

reader = ReaderRay.RayReader(OUTPUT_PATH, READ_VALUE, HORIZONTAL, HORIZONTAL_GRID, SERIES, SERIES_GRID, MOVIE, MOVIE_GRID)

reader.set_analytical(ANALYTICAL)
reader.set_auxiliary(AUXILIARY)
reader.read()

reader.set_ylimits(LIMITS)
reader.set_delay(0.1)

reader.plot(numerical=True, analytical=True, auxiliary=False)