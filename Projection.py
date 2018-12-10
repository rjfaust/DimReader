from DualNum import*
paramOpts = []
class Projection:
    def __init__(self, points, parameters):
        self.origPoints = points
        self.projectionParams = parameters

    def run(self):
        return self.origPoints

    def getExecParams(self):
        return None

    def setExecParams(self,params):
        pass
