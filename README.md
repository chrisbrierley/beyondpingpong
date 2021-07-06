# beyondpingpong
A possible thing to look at during the Heidelberg Workshop.

The science idea is that we can assess whether the power within a spectral band is correct in a model. Can we then pull that knowledge back into the (non-forced) component of a future projection? Here I've used a single model's with initial condition ensembles (this simplifies the situation by removing the need to consider model uncertianty). I have little idea why the past multi-decadal climate in the Last Millennium Reanalysis is less variable than that simulated by CESM. It seems to be counter to the suggestion that models are generally insufficiently variable, which would require to greater spread in future predictions.   

The codes and data in this repository will result in something like the figure below. The NCL code was used to do the actual analysis, then Inkscape was just to make the plot look cleaner. I've tried to translate it to a python notebook, but it hasn't got anywhere near as far and doesn't yet have the required bandpass filtering included. 


![Possible Image](Pretty_image.png)

*Pilot work demonstrating the potential for projection variance calibration using the past 1000 years. (a) Annual, global mean temperature anomalies in from the last millennium reanalysis and simulations. (b) Distribution of multidecadal anomalies in both datasets. (c) Multidecadal projections before and after the variance calibration.*
