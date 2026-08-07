import matplotlib.pyplot as plt

class Data3D:

    def __init__(self, var_name: str, x_var: str, x_grid, y_var: str, y_grid, z_var: str, z_grid, c):
        self.var_name = var_name
        self.x_var = x_var.lower()
        self.x_grid = x_grid
        self.y_var = y_var.lower()
        self.y_grid = y_grid
        self.z_var = z_var.lower()
        self.z_grid = z_grid
        self.c = c

        self.xlimits = None
        self.ylimits = None
        self.legend_prefix = ''

        self.x_val = []
        self.y_val = []
        self.z_val = []

        self.num_plot = [[[]]]
        self.an_plot = [[[]]]
        self.aux_plot = [[[]]]

        self._delay = 0.01


    def set_xlimits(self, limits):
        self.xlimits = limits


    def set_ylimits(self, limits):
        self.ylimits = limits
    

    def set_legend_prefix(self, prefix):
        self.legend_prefix = prefix


    def set_delay(self, delay):
        self._delay = delay
    

    def integrate_plot(self, plot, d_ax='x', start=0):
        if d_ax.lower() == 'x':
            for z_iter in range(len(self.z_grid)):
                for y_iter in range(len(self.y_grid)):
                    for x_iter in  range(start+1, len(self.x_grid)):
                            del_x = self.x_val[x_iter] - self.x_val[x_iter - 1]
                            plot[z_iter][y_iter][x_iter] = plot[z_iter][y_iter][x_iter - 1] + plot[z_iter][y_iter][x_iter] * del_x
        elif d_ax.lower() == 'y':
            for z_iter in range(len(self.z_grid)):
                for y_iter in range(start+1, len(self.y_grid)):
                    for x_iter in  range(len(self.x_grid)):
                            del_y= self.y_val[y_iter] - self.y_val[y_iter - 1]
                            plot[z_iter][y_iter][x_iter] = plot[z_iter][y_iter - 1][x_iter] + plot[z_iter][y_iter][x_iter] * del_y
        elif d_ax.lower() == 'z':
            for z_iter in range(start+1, len(self.z_grid)):
                for y_iter in range(len(self.y_grid)):
                    for x_iter in  range(len(self.x_grid)):
                            del_z = self.z_val[z_iter] - self.z_val[z_iter - 1]
                            plot[z_iter][y_iter][x_iter] = plot[z_iter - 1][y_iter][x_iter] + plot[z_iter][y_iter][x_iter] * del_z
        else:
            print(f"Intrgration axis {d_ax.lower()} not recognized as 'x', 'y', or 'z'")


    def abs_integrate_plot(self, plot, d_ax='x', start=0):
        if d_ax.lower() == 'x':
            for z_iter in range(len(self.z_grid)):
                for y_iter in range(len(self.y_grid)):
                    plot[z_iter][y_iter][0] = abs(plot[z_iter][y_iter][0])
                    for x_iter in  range(start+1, len(self.x_grid)):
                            del_x = self.x_val[x_iter] - self.x_val[x_iter - 1]
                            plot[z_iter][y_iter][x_iter] = abs(plot[z_iter][y_iter][x_iter])
                            plot[z_iter][y_iter][x_iter] = plot[z_iter][y_iter][x_iter - 1] + plot[z_iter][y_iter][x_iter] * del_x
        elif d_ax.lower() == 'y':
            for z_iter in range(len(self.z_grid)):
                for x_iter in range(len(self.x_grid)):
                    plot[z_iter][0][x_iter] = abs(plot[z_iter][0][x_iter])
                    for y_iter in  range(start+1, len(self.y_grid)):
                            plot[z_iter][y_iter][x_iter] = abs(plot[z_iter][y_iter][x_iter])
                            del_y= self.y_val[y_iter] - self.y_val[y_iter - 1]
                            plot[z_iter][y_iter][x_iter] = plot[z_iter][y_iter - 1][x_iter] + plot[z_iter][y_iter][x_iter] * del_y
        elif d_ax.lower() == 'z':
            for x_iter in range(len(self.x_grid)):
                for y_iter in range(len(self.y_grid)):
                    plot[0][y_iter][x_iter] = abs(plot[0][y_iter][x_iter])
                    for z_iter in  range(start+1, len(self.z_grid)):
                            plot[z_iter][y_iter][x_iter] = abs(plot[z_iter][y_iter][x_iter])
                            del_z = self.z_val[z_iter] - self.z_val[z_iter - 1]
                            plot[z_iter][y_iter][x_iter] = plot[z_iter - 1][y_iter][x_iter] + plot[z_iter][y_iter][x_iter] * del_z
        else:
            print(f"Intrgration axis {d_ax.lower()} not recognized as 'x', 'y', or 'z'")
    

    def plot_z(self, z, ax, numerical, analytical, auxiliary):                    
        """Plots information at a set z value given an pyplot axis
        Does not actually display the plot
        """

        for y_iter in range(len(self.y_grid)):
            if len(self.y_grid) > 1:
                label = f'{self.legend_prefix}{self.y_var}={round(self.y_grid[y_iter], 2)}'
            else:
                label = self.legend_prefix

            if numerical:
                ax.scatter(self.x_val, self.num_plot[z][y_iter], label=f'num {label}', marker='.', c=self.c)

            if analytical:
                ax.plot(self.x_val, self.an_plot[z][y_iter], label=f'an {label}', c=self.c)
                    
            if auxiliary:
                ax.plot(self.x_val, self.aux_plot[z][y_iter], label=f'aux {label}', marker='x', c=self.c)

        ax.set(xlabel=f'{self.x_var}', ylabel=f'{self.var_name}')
        plt.title(f"{self.z_var} = {round(self.z_val[z], 2)}")
    
    def plot(self, numerical=True, analytical=False, auxiliary=False):
        "Runs animated plots based off the values found"
        _, ax = plt.subplots(1)
        for z_iter in range(len(self.z_grid)): #type: ignore
            ax.clear()
            self.plot_z(z_iter, ax, numerical, analytical, auxiliary)
            ax.set_xlim(self.xlimits)
            ax.set_ylim(self.ylimits)
            plt.legend()
            plt.pause(self._delay)
        plt.show()