CoxHD
-----
An R package with high-dimensional extensions of the Cox proportional hazards model.

It includes:

####CoxRFX: A random effects model
This model assumes a shared normal distribution of effect sizes. The data can be partitioned into groups with different effect mean and variance. 
It can be useful for fitting $n\approx p$ dimensions.
 
####CoxCPSS: Complementary pairs stability selection
Stability selection is used for variable selection. Error control is possible based on the log-concave model by Shah and Samworth.

##Installation

Installation is easy using devtools::install_github()

	> library(devtools); install_github("mg14/CoxHD/CoxHD")