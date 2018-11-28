### Software Requirements ###

* Install sage:
  * 8.2 version (released on 2018-05-05) is needed. (http://www.sagemath.org/)

* Add SageMath directory to $PATH

* Install ore_algebra: follow the instructions in http://www.kauers.de/software.html

### How to Run Codes ###

* go to 'code/examples' folder.

type the command line
```
~/Downloads/SageMath/sage exampleForKrawczykTest.sage
```
to execute the scripts

* or in 'code/examples' folder, run sage by the command line
```
~/Downloads/SageMath/sage
```

When the sage is running, load the script. For example
```
load("expExample")
```


* Scripts in 'code/examples' folder.
     * err-2equations-comparison.sage : compare alpha theory, Krawczyk method based tests on the system with 2 error functions.(example in section 5.1)
     * mixed-2equations.sage : computes gamma value using maximum value of a Bessel function and an error function.(example in section 5.2)
     * expExample.sage : comparison between 'alphaCertified' on the system with an exponential function.(example in section 5.3)
     * ellipseExp-alpha.sage : maximization problem with two ellipses certified by alpha theory (example in section 5.4)
     * ellipseExp-Krawczyk.sage : maximization problem with two ellipses certified by Krawczyk mehtod (example in section 5.4)

