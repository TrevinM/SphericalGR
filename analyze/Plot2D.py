import matplotlib.pyplot as plt
import ReaderCrit


PATH = '../test/CRIT_ODEq/iteration_9'

X_NAME = 'time'
X_COL = 0

Y_NAME = 'min lapse'
Y_COL = 6

reader = ReaderCrit.CritReader(PATH, '.mon', X_NAME, X_COL, Y_NAME, 7)

reader.read()
reader.set_ylimits([0,1])
reader.set_legend_prefix('eta')
reader.plot()

reader2 = ReaderCrit.CritReader(PATH, '.mon', X_NAME, X_COL, Y_NAME, 6)
reader2.read()
reader2.set_ylimits([0,1])
reader2.set_legend_prefix('eta')
reader2.plot()

plt.legend()
plt.show()