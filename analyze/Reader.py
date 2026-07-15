import glob
from ThreeDPlotter import ThreeDPlotter
from abc import ABC, abstractmethod

class Reader(ThreeDPlotter, ABC):

    def __init__(self, var_name: str, x_var: str, x_grid, y_var: str, y_grid, z_var: str, z_grid):
        super(Reader, self).__init__(var_name, x_var, x_grid, y_var, y_grid, z_var, z_grid)

        #Abstract
        self._files: list

        #Protected
        self._analytical = None
        self._aux = None

        #Public
        self.resolution = self._grid_resolution()
    
    
    @abstractmethod
    def _grid_resolution(self):
        """Reads the name of the first file found to find information about grid resolution
        Returns: [n_r, n_theta, n_time]
        """
        raise NotImplementedError


    @abstractmethod
    def _find_time(self, lines):
        """Given a file read using readlines will return the current time value"""
        raise NotImplementedError
    

    @abstractmethod
    def _find_value(self, lines, r_step, th_step):
        """Given a file read using readlines and r_step, th_step returns dumped values
        Returns: [r_val, th_val, func_val]
        """
        raise NotImplementedError

    
    def set_analytical(self, func):
        """Define an analytical function that can be evaluated along with the numerical slices
        func : A function f(r, theta, time)
        """
        self._analytical = func


    def set_auxiliary(self, func):
        """Define an auxiliary function that can compare the numerical and analytical
        func : A function f(r, theta, time)
        """
        self._aux = func


    def read(self):
        "Populate the reader with variable information"
        self._initialize_axes()
        self._initialize_data()


    def _initialize_axes(self):
        self.r_grid = []
        self.theta_grid = []
        self.time_grid = []

        self.x_index = 0

        if self.x_var == 'r':
            self.x_index = 0
            if self.x_grid == None:
                self.x_grid = range(self.resolution[self.x_index])
            self.r_grid = self.x_grid
        elif self.x_var == 'theta':
            self.x_index = 1
            if self.x_grid == None:
                self.x_grid = range(self.resolution[self.x_index])
            self.theta_grid = self.x_grid
        elif self.x_var == 't':
            self.x_index = 2
            if self.x_grid == None:
                self.x_grid = range(self.resolution[self.x_index])
            self.time_grid = self.x_grid
        else:
            print(f"x variable {self.x_var} not recognized as 'r', 'theta' or 't'.")

        print(f"x index: {self.x_index}")

        self.y_index = 0

        if self.y_var == 'r':
            self.y_index = 0
            self.r_grid = self.y_grid
        elif self.y_var == 'theta':
            self.y_index = 1
            self.theta_grid = self.y_grid
        elif self.y_var.lower() == 't':
            self.y_index = 2
            self.time_grid = self.y_grid
        else:
            print(f"y variable {self.y_var} not recognized as 'r', 'theta' or 't'.")

        print(f"y index: {self.y_index}")

        self.z_index = 0

        if self.z_var == 'r':
            self.z_index = 0
            if self.z_grid == None:
                self.z_grid = range(self.resolution[self.z_index])
            self.r_grid = self.z_grid
        elif self.z_var.lower() == 'theta':
            self.z_index = 1
            if self.z_grid == None:
                self.z_grid = range(self.resolution[self.z_index])
            self.theta_grid = self.z_grid
        elif self.z_var == 't':
            self.z_index = 2
            if self.z_grid == None:
                self.z_grid = range(self.resolution[self.z_index])
            self.time_grid = self.z_grid
        else:
            print(f"z variable {self.z_var} not recognized as 'r', 'theta' or 't'.")

        print(f"z index: {self.z_index}")

        print(f"time grid: {self.time_grid}")
        print(f"r grid: {self.r_grid}")
        print(f"theta grid: {self.theta_grid}")


    def _initialize_data(self):
        self.num_plot = [[[] for _ in range(len(self.y_grid))] for _ in range(len(self.z_grid))] # type: ignore
        self.an_plot = [[[] for _ in range(len(self.y_grid))] for _ in range(len(self.z_grid))] # type: ignore
        self.aux_plot = [[[] for _ in range(len(self.y_grid))] for _ in range(len(self.z_grid))] # type: ignore

        for k, t_step in enumerate(self.time_grid):
            f = open(self._files[t_step])
            print(f"Reading {self._files[t_step]}")
            lines = f.readlines()

            t = self._find_time(lines)
            
            for i, r_step in enumerate(self.r_grid):
                for j, th_step in enumerate(self.theta_grid):
                    ijk = [i, j, k]
                    x_iter = ijk[self.x_index]
                    y_iteration = ijk[self.y_index]
                    z_iteration = ijk[self.z_index]

                    vals = self._find_value(lines, r_step, th_step)
                    r = vals[0]
                    th = vals[1]
                    num_val = vals[2]

                    self.num_plot[z_iteration][y_iteration].append(num_val)

                    #Run analytical and auxiliary functions if given
                    if self._analytical != None:
                        an_val = self._analytical(r, th, t) #type: ignore
                        self.an_plot[z_iteration][y_iteration].append(an_val)   

                        if self._aux != None:
                            aux_val = self._aux(num_val, an_val) #type: ignore
                            self.aux_plot[z_iteration][y_iteration].append(aux_val)

                    # Give values to the axes
                    if y_iteration == 0:
                        if z_iteration == 0:
                            if (self.x_index == 2):
                                self.x_val.append(t)
                            else:
                                self.x_val.append(vals[self.x_index])
                        if x_iter == 0:
                            if (self.z_index == 2):
                                self.z_val.append(t)
                            else:
                                self.z_val.append(vals[self.z_index])  
                                print(f"z_val appended {self.z_val[-1]}")
