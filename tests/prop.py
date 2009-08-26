class C(object):

    def __init__(self):
        self.__x = 0

    def getx(self):
        return self.__x

    def setx(self, x):
        if x < 0: x = 0
        self.__x = x

    x = property(getx, setx, doc = "I am doc")

a = C()

a.x = 10
print a.x

print a.x.__doc__
