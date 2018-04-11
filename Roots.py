import math

print('\n',"Numerical methods for finding roots of nonlinear algebraic and transcendental equations",'\n')

A = -5.0; B = 5.0; eps = 6; e = 10**(-eps); ans = 15
def f(x):
    return 1.2*(x**4)+2*(x**3)-13*(x**2)-14.2*x-24.1
def fp(x):
    return 4.8*(x**3)+6*(x**2)-26*x-14.2

#RootSep
n = 50.0
h = (B-A)/n
x1 = A; x2 = x1+h
k = 0
ints = []
while(x2 <= B):
    y1 = f(x1)
    y2 = f(x2)
    if (y1*y2 <= 0):
        k += 1
        #ints.append((x1,x2))
        ints.append((round(x1,eps),round(x2,eps)))
    x1 = x2
    x2 += h
print('\n',"Number of intervals:", k)
for i in range(0,k):
    print('[',ints[i][0],',',ints[i][1],']')
print('\n')

#We got
print("[A,B] = [",A,',',B,']','\n',"h = ",h,'\n',"eps = ",e,'\n',"f(x) = 1.2x^4+2x^3-13x^2-14.2x-24.1",'\n',sep="")

#Bisection
print("----------------------",'\n',"Bisection method")
for i in range(0,k):
    a = ints[i][0]; b = ints[i][1]
    x = (a + b) / 2
    print("For x", i, ' initial approximation is: ', x, sep="")
    count = 1
    while math.fabs(f(x)) >= e:
        a, b = (a, x) if f(a) * f(x) < 0 else (x, b)
        x = (a + b) / 2
        count += 1
    print("x",i,' = ',x, sep="")
    print("|f(x", i, ')-0| = ', (math.fabs(f(x))),'\n', sep="")
    print('\n',"With",count,"iterations.",'\n')

#Newton
print("----------------------",'\n',"Newtons method")
for i in range(0,k):
    a = ints[i][0]; b = ints[i][1]
    x0 = (a + b) / 2
    print("For x", i, ' initial approximation is: ', x0, sep="")
    x1 = x0 - (f(x0) / fp(x0))
    count = 1
    while (math.fabs(x1 - x0) >= e):
        x0 = x1
        x1 = x0 - (f(x0) / fp(x0))
        count += 1
    print("x",i,' = ',x1, sep="")
    print("|f(x", i, ')-0| = ', round(math.fabs(f(x1)),ans),'\n', sep="")
    print('\n',"With",count,"iterations.",'\n')

#Modified Newton
print("----------------------",'\n',"Modified Newtons method")
for i in range(0,k):
    a = ints[i][0]; b = ints[i][1]
    x0 = (a + b) / 2
    print("For x", i, ' initial approximation is: ', x0, sep="")
    fp0 = fp(x0)
    x1 = x0 - (f(x0) / fp0)
    count = 1
    while (math.fabs(x1 - x0) >= e):
        x0 = x1
        x1 = x0 - (f(x0) / fp0)
        count += 1
    print("x",i,' = ',x1, sep="")
    print("|f(x", i, ')-0| = ', round(math.fabs(f(x1)),ans),'\n', sep="")
    print('\n',"With",count,"iterations.",'\n')

#Secant
print("----------------------",'\n',"Secant method")
for i in range(0,k):
    a = ints[i][0]; b = ints[i][1]
    x0 = a; x1 = b
    x = x0 - f(x0) * (x1 - x0) / (f(x1) - f(x0))
    print("For x", i, ' initial approximation is: ', x, sep="")
    count = 1
    while (math.fabs(x - x1) >= e):
        x0 = x1; x1 = x
        x = x0 - f(x0) * (x1 - x0) / (f(x1) - f(x0))
        count += 1
    print("x", i, ' = ', x, sep="")
    print("|f(x", i, ')-0| = ', round(math.fabs(f(x)),ans),'\n', sep="")
    print('\n',"With",count,"iterations.",'\n')

#Iteration
print("----------------------",'\n',"Iteration method")
for i in range(0,k):
    a = ints[i][0]; b = ints[i][1]
    x1 = b; lam = 1 / fp(a)
    x = x1 - lam * f(x1)
    print("For x", i, ' initial approximation is: ', x, sep="")
    count = 1
    while (math.fabs(x - x1) >= e):
        x0 = x1; x1 = x
        x = x1 - lam * f(x1)
        count += 1
    print("x", i, ' = ', x, sep="")
    print("|f(x", i, ')-0| = ', round(math.fabs(f(x)),ans),'\n', sep="")
    print("With",count,"iterations.",'\n')
