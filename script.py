for i in range(0, 3):
    print '''
torus 0 0 0 10 400 0.5
torus 0 0 0 15 400 0.4
torus 0 0 0 20 400 0.3
torus 0 0 0 25 400 0.2
torus 0 0 0 30 400 0.1
torus 0 0 0 35 400 0.05
torus 0 0 0 40 400 0.01
x 0
y %d
a
i
sphere 0 0 0 50 0.5
sphere 0 0 0 75 0.4
sphere 0 0 0 100 0.3
sphere 0 0 0 125 0.2
sphere 0 0 0 150 0.1
sphere 0 0 0 200 0.05
sphere 0 0 0 250 0.01
x 30
y 30
a
i
f imgs/frame%d.png
clear''' % (i*3,i)
