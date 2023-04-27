Mapping cells through time, space, and beyond
=============================================
 <img src="images/light_mode_logo.png"
   alt="Markdown Monster icon"
   style="float: left; margin-right: 10px;" />

   [`moscot`](moscot-tools.org) is a toolbox for Optimal Transport (OT) applications in single-cell genomics. Moscot makes OT applications accessible to biologists by providing a unified framwork. Being built upon [OTT-JAX](https://ott-jax.readthedocs.io/en/latest/index.html), moscot leverages state-of-the-art advances in OT algorithms allowing for the application to large-scale single-cell datasets while enabling the incorporation of multiple modalities. While the documentation of the [software package](https://github.com/theislab/moscot) can be found on [moscot-tools.org](moscot-tools.org), this repository contains the data analyses perfomed for the _preprint.

Manuscript
----------
Our manuscript is available as a [preprint](biorxiv.org) on bioRxiv. 

Reproducibility
---------------
To ease reproducibility of our preprint results, the structure of this repository follows the one of the manuscript.
Each folder contains notebooks and scripts necessary to reproduce the corresponding analysis. 
Note that benchmarks are in another [repository](https://github.com/theislab/moscot_benchmarks).

Results
-------
The analysis can be found in the following directories

 - moscot.time (Fig. 2): `analysis/time`
 - moscot.space (Fig. 3): `analysis/space`
 - moscot.spatiotemporal (Fig. 4): `analysis/spatiotemporal`
 - pancreas multiome (Fig. 5): `analysis/pancreas`

