# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from Node import *
from Material import *
from Element import *
#from Truss import *
#from FiberMaterial import *
#from System import *
def Problem1():
    # Use a breakpoint in the code line below to debug your script.
    a = Node(0,0,0)
    b = Node(0,0,1)
    mat = Material()
    elem = Element(a,b,mat )
    print(a)
    print(a.__repr__())# Press Ctrl+F8 to toggle the breakpoint.
    print(mat)
    print(mat.__repr__())  # Press Ctrl+F8 to toggle the breakpoint.
    print(elem)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    Problem1()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
