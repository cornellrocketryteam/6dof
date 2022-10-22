#Unscented Kalman Filter Implimentation for Cornell Rocketry Real-Time State Estimation
#Developed by Sam Noles in 2022

#PyPlot Plotting Examples: https://gist.github.com/gizmaa/7214002#multiaxis

#########info##############

#system state: z = [r v q w]
#r: position of R (tip of nose cone) wrt to O (ground) in inertial frame 1,2,3
#v: inertial frame derivative of r in intertial frame coordinates 4,5,6
#q: quaternion representing rotation from inertial to body frame 7,8,9,10
    #aka apply q to a vector written in I cordinates to get the same vector written in body frame coordinates
#w:IwB rotatoin of B frame wrt to the inertial frame 11,12,13

include("6dof.jl")



#functions 

