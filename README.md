# IRLS plane fitting
Iterative Reweighted Least-squares (IRLS) using M-estimator for plane fitting 

### Description
The main idea behind this project is to have a robust plane fitting algorithm that will be less susceptible to noise compared to a standard plane fitting approach, for example `ax + by + cz = d`. While the standard least-squares method tries to minimize sum of the distances between points to the fitted plane, IRLS tries to minimize weighted sum of the distances, thus reducing effect of outliers.

### Current Functionality
Build and run `./irls_plane_fitting`

### Example Results
![image1](https://github.com/ysimonov/IRLS_plane_fitting/blob/main/example/result_1.png)
![image2](https://github.com/ysimonov/IRLS_plane_fitting/blob/main/example/result_2.png)
![image3](https://github.com/ysimonov/IRLS_plane_fitting/blob/main/example/result_3.png)
