In this study, a real-time model updating method is developed that uses the local eigenvalue modification procedure (LEMP) to break down the original eigenvalue solution into a collection of second-order equations. The suggested technique was developed offline after experimental validation utilizing data from the DROPBEAR testbed.

# Model Development

![flowchart](https://user-images.githubusercontent.com/69466658/183703072-97f72333-7a72-4f6f-bc6a-352f96b68397.PNG)

The objective of the experimental process is to ascertain the actual system response of the DROPBEAR testbed under various boundary circumstances. This is done by using the accelerometer that is positioned on the free end of the beam to gather acceleration data from the system. A sliding Hann window is used to smooth the time-series data before the acceleration data are analyzed. The Fast Fourier Transform (FFT) of the acceleration data is then used to determine the beam's inherent frequency. Then, using comparison methods, the observed system response is compared to the several analytically calculated models to estimate the state. The analytical step is to create model in which LEMP is implemented on selected roller locations to calculate the natural frequencies of the beam. The calculate frequencies are then compared to the measured frequecies using two comparison criteria, error minimization and bounded regression.

# Model updating results
Error minimization and linear regression were used as comparison techniques in the calculation of roller estimations by GE solutions, LEMP solutions, and LEMP solutions in a Bayesian search space.

![error_minimization2](https://user-images.githubusercontent.com/69466658/183679692-5af63bb4-afc1-4b92-80cb-8a89e4e9ab6f.jpg)

for the figure above, when using the error minimization strategy, the roller would move and fluctuate, producing abrupt transitory errors in most cases. But the estimations improved with the use of a Bayesian search space. When using the bounded regression technique, estimation fluctuation decreased during roller movement but not during stationary times.

As the number of nodes used in the model increases, the estimation became better and less error is observed in the model update results as shown below.

![LEMP_error](https://user-images.githubusercontent.com/69466658/183681404-0d3e561f-1113-41f0-8bb6-a2c6c89c6ca9.JPG)

# Model timing and accuracy

*Solver time*

![solver_time](https://user-images.githubusercontent.com/69466658/183702159-6388887f-ffb8-4836-b061-ba475306598b.jpg)

![update_time](https://user-images.githubusercontent.com/69466658/183702198-44da4b6d-04ce-46b4-bc68-78cabc0e6235.jpg)

*Accuracy and optimal configuration*

![results_table](https://user-images.githubusercontent.com/69466658/183683777-b6fd04fb-3eaf-4929-ad8d-29edf04ff306.PNG)

![MEA vs iteration time](https://user-images.githubusercontent.com/69466658/183682900-38fbadb6-edfe-4215-bfd6-8a5df3d0ea84.jpg)
