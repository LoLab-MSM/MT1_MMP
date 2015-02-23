import numpy
a = numpy.array([1, 0, 2])
a2 = numpy.array([10, 11, 12])
a1 = [numpy.array([2,7,9]),numpy.array([6,1,0]),numpy.array([0,7,4])]
b1 = numpy.array([numpy.array([1,0,0]),numpy.array([0,1,0]),numpy.array([0,0,1])])

t = a1[0] * a
a3= numpy.array([numpy.array([1,2,3]), numpy.array([4,5,6]), numpy.array([7,8,9])])
rr = a3 * b1
rr2 = b1 * a3
r = a*a2
b = [[3], [0]] + [[0], [2]]
c = a * b
d = a*[[3], [0], [1]] + a*[[0], [2], [9]]
a1 = [numpy.array([2,7,9]),numpy.array([6,1,0]),numpy.array([0,7,4])]
M=numpy.identity(3)
print M[0]
G = a1[0] * M[0] + a1[1] * M[1]+ a1[2] * M[2] 
print G
# 
# print c
print rr
print rr2
# print M
# print r
# print a3
# print t