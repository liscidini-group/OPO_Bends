
# Bends

This notebook contains the formulas to obtain the coordinates of different kinds of curves. Here you can find the comments for each function and how to use it.
In the Jupyter Notebook <code>OPO_curves</code>  you can find the functions and the plots of the curves we used for the design of Linearly Uncoupled Resonator for OPO. 


## Overview example

In the Jupyter Notebook <code>Euler'sbend&co.ipynb</code>  you can find a brief explanation of the used formulas and the plots for each curve (considering generic parameters). 
In the notebook <code>Sbend.ipynb</code> you will find a more detailed description of the S-bends, with comparison between different types of S-bends, along with the functions to plot a series of S-bends. You can 

With this code you can obtain the coordinates of the following bends:

 1. U-bends 
- Non simmetrical
- Total Euler
- Hybrid Euler-circular
- 4-Euler 180
- Cubic Bezier

 2. S-bends
 - Sinusoidal
 - Double Euler
 - Four Euler
 - Cubic Bezier

 3. 90° bends
- Euler
- N-adjustable Euler
- Cubic Bezier



You can use all the modules by adding the following line at the beginning of your script, where you import other libraries:

```py
import bend_functions as bf
```

Each functions returns two vectors of coordinates (x,y) for each curves.


# U-bends

U bends can be realized in different ways. 

### Non symmetrical U-bend

The easiest way is to consider the expression for the Euler bend or N-adjustable Euler bend and fix the final angle of the curve to be $\pi$. 
The functions to be used are

```py
euler_bend(R, theta, num_points)
```

```py
euler_bend_NA(R, theta, N, num_points)
```

where R is the final radius of curvature for the Euler, $\theta$ is the total angle ($\theta=\pi$ in this case), N is the order of the NA bends (for N=1 we re-obtain the Euler bend) and num_points is the number of coordinates of the final vectors.
These are not symmetric curves and lead to longher path distance.

### Hybrid Euler-circular U-bends
To obtain a simmetric U-bend we interconnect two Euler's bend with an arc of circumference, as in [Ji, Xinru, et al.,Communications Physics 5.1 (2022): 84.].
The function used is 
```py
euler_bend_Hybrid_180(Reff, p, num_points):
```

where $R_{eff}$ is the radius of the relative circular curve (or half the distance of the initial and final point of the curve). p is a parameter such that $0\leq p\leq 1$, (p=0 is the circular curve, p=1 is the total Euler curve). num_points fixed the number of points of the curves.

Another type of 180 degree bend is also considered. It is formed by 4 Euler, from 0 to $\pi$/4, from $\pi/4$ to $\pi/2$ and so on. Basically they correspont to a 90 degrees symmetric Euler. The function is
```py
euler_bend_NA_180(Reff, n, num_points):
```

where R is the effective radius and n is the order of the Euler NA bend.

### Bezier U-bends
Finally we consider also the expression for the Cubic bezier. The Bezier curve are defined through a class and thay can be obtain with the functions
```py
curve = bf.CubicBezier(p0x=0, p0y=0, p1x=2*R-B2, p1y=B1, p2x=2*R-B2, p2y=2*R-B1, p3x=0, p3y=2*R)
x,y= curve.calc_curve()
```

The coordianates of the 4 points have to be inserted. They can be written by fixing the first and the last point to (0,0) and (0,$2R_{eff}$), and changing the two central points in function of the parameters B1 and B2.



# S-bends

Let us now consider different expression for the S-bends. We call d the lateral displacement of the curve. Each curve has a parameter that can be varied to adjust the curvature of the bend. Each curve is traslated such that the beginning is fixed at (0,0).


### Sinusoidal functions
The considered equation and the respective functions are

$$ y= \dfrac{d}{2}\tanh{\left(\dfrac{x}{T}\right)}$$

```py
tanh(T,d,L,npts)
```


$$ y= \dfrac{d}{2}\mathrm{erf}(x/E) $$
```py
erf(E,d,L,npts)
```


$$y=\dfrac{2d}{\pi}\arctan(\pi C x/2)$$
```py
arctan(C,d,L,npts)
```


$$y=\dfrac{dCx^2}{2\sqrt{1+(Cx)^2}}$$
```py
rad(C,d,L,npts)
```


$$ y=\dfrac{dxC}{2(1+|xC|)} $$
```py
abv(C,d,L,npts)
```

Where d is the lateral displacement, L is the forward length and each function has a parameter (T,E,C,...) that can be varied to adjust the curvature of the S-bend. Each function has initial points set to (0,0) and final to (L,d)

### Euler S-bend

With the function
```py
euler_bend_S2(theta, R, d,  n,num_points)
```
is possible to obtain the coordinates for the S-bend obtained by the connection of 2 Euler's bend. The curvature at the central point can be controlled by the parameter $\theta$. R is the minimum radius of curvature, d is the lateral displacement and n is the order of the N-adjustable Euler's bend.

The function 
With the function
```py
euler_bend_S4(theta, R, d,  n,num_points):
```
allows to obtain a S-bend made of 4 Euler. In this way, the discontinuity at the center of the bend is eliminated. The parameter $\theta$ can be controlled to adjust the central slope. N control the order of the N-adjustable Euler bend. R, d and n hace the same meaning for the 2-Euler function.

### Bezier S-bend

With the functions
```py
curve = bf.CubicBezier(p0x=0, p0y=0, p1x=L/2-B, p1y=0, p2x=L/2+B, p2y=d, p3x=L, p3y=d)
x,y= curve.calc_curve()
```
the Bezier bend can be obtained and ajdusted by controlling the parameter B, for different values of d and L.

# 90° bends

### Euler and NA 90° bends
The used function is 
```py
euler_bend_NA_90(Reff, n, num_points):
```
Where Reff fixes the effective radius of the bend, n is the order of the N-adjustable Euler's bend and num_points is the number of the coordinates for each curve.

### Bezier 90° bends 
The used functions are
```py
curve = bf.CubicBezier(p0x=0, p0y=0, p1x=(1-B)*R, p1y=0, p2x=R, p2y=B, p3x=R, p3y=R)
x,y= curve.calc_curve()
```

where the coordinates of the two central points are controlled by varying the parameter B.



## Big rings
In the  Jupyter Notebook <code>Big_rings.ipynb</code>  you can find a specific application of these functions to the design of two big rings.
