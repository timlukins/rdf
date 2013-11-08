Registered Differential Features
=====

"A discretised set of 15 categorical labels to describe the changes in curvature over time."[1]

***

Table of Contents
-----------------

**[Introduction](#introduction)**

**[Theory](#theory)**

**[Implementation](#implementation)**

**[Application](#application)**

**[References](#references)**

***

Introduction
------------

Imagine a thin sheet of metal. Think of all the possible ways it could be distorted without tearing it... Then look at this diagram:   

![Graph of 15 curvature changes](graph.png) 

In essense, there are *15 possible ways* that the surface can distort::

** Those that typify, and exemplify further, the formation of a prototype form (e.g. *protrude* a peak, *subside* a pit, *fold* a valley, *bend* a ridge, *warp* a saddle).

** Those that move the opposite way from the prototype form towards flat (e.g. *flatten*). In combination with the first transition, this serves to capture the equivalent variation in curvedness expressed by Koenderink.

** Those bi-directional transitions that are also applicable only between neighbouring non-flat prototypes (e.g. to *squeeze* a pit to form a valley, to *collapse* a valley to form a pit, to *dimple* a valley to form a saddle, to *crumple* a saddle to form a valley, to *crease* a saddle to form a ridge, to *dent* a ridge to form a saddle, to *bulge* a ridge to form a peak, and *stretch* a peak to form a ridge).
  
** Those shapes that do not make any transition - as they have no observable change in curvature (e.g. they are *constant*).

This results in a total of *15 different deformation classes*, which we formalise as the *type* of deformation.


Theory
------

<!-- Use http://www.url-encode-decode.com/urlencode to encode equations -->

More precisely, consider the two principal curvatures centred on a surface neighbourhood. The principal curvatures are the two maximal changes in the form of the surface around that point. They have a direction and a magnitude of change. And they are orthogonal - i.e. at right angles to one another.

Now, consider how these two properties - the respective magnitutes of the principle curvatures - change w.r.t. one another. 

On this basis, it can be realised that there are only a limited number of different types of transition which we enumerate in more detail here:

** Those that typify, and exemplify further, the formation of a prototype (e.g. *protrude* a peak, *subside* a pit, *fold* a valley, *bend* a ridge, *warp* a saddle).

** Those that move the opposite way from the prototype towards flat (e.g. *flatten*). In combination with the first transition, this serves to capture the equivalent variation in curvedness expressed by Koenderink.

** Those bi-directional transitions that are also applicable only between neighbouring non-flat prototypes (e.g. to *squeeze* a pit to form a valley, to *collapse* a valley to form a pit, to *dimple* a valley to form a saddle, to *crumple* a saddle to form a valley, to *crease* a saddle to form a ridge, to *dent* a ridge to form a saddle, to *bulge* a ridge to form a peak, and *stretch* a peak to form a ridge).
  
** Those shapes that do not make any transition - as they have no observable change in curvature (e.g. they are *constant*).

This results in a total of *15 different deformation classes*, which we formalise as the *type* of deformation:

    T\in[1,...,15]\label{eq:type}

that can occur over any given duration as defined by {\it the relative change in the principal curvatures $\Delta\kappa_{1}$ and $\Delta\kappa_{2}$.

As with the initial shape classes, in order to define the zero boundary region it is necessary to employ a threshold term ($\Delta\kappa=0\iff-\theta_{change}<\Delta\kappa<\theta_{change}$). 

We furthermore can define the *extent of change $(E)$ to also measure the degree to which this deformation occurs over the duration. 

    E=\sqrt{\frac{\Delta\kappa_{1}^{2}+\Delta\kappa_{2}^{2}}{2}}.

![equation](http://latex.codecogs.com/gif.latex?1%2Bsin%28mc%5E2%29%0D%0A)

(denoted $\kappa_1$ and $\kappa_2$).

S=\frac{2}{\pi}arctan\left((\kappa_{1}+\kappa_{2})/(\kappa_{1}-\kappa_{2})\right)

C=\sqrt{(\kappa_{1}^{2}+\kappa_{2}^{2})/2}.

Implementation
-----------

```python
import cv 
import cv2
import numpy as np 

```

Application
-----------




References
----------

http://homepages.inf.ed.ac.uk/rbf/PAPERS/lukins06qualitative.pdf

A. M. McIvor and R. J. Valkenburg. Principal frame and principal quadric estimation. Image and Vision Computing New Zealand, pages 55–60, 1996

J. J. Koenderink and A. J. van Doorn. Surface shape and curvature scales. Image and Vision Computing, 10(8):557–
565, 1992.

http://www.eecs.berkeley.edu/~jima/mypapers/AS_CADA2013_quadrics.pdf
