from Projection import*

class TangentMapProjection(Projection):

    def __init__(self, points, parameters):
        self.origPoints = points
        params = [self.filename]

        for i in range(len(parameters)):
            params[i] = parameters[i]