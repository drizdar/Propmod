class membrane:
    def __init__(self, dimensions, properties):
        self.dims = dimensions
        self.props = properties
    def discretise(self, n):
        self.n_evaluation_points = n
        
