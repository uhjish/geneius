

class GeneiusError:
    
    def __init__(self,msg):
        self.message = msg
    
    def __repr__(self):
        return "GeneiusError: %s" % self.message
