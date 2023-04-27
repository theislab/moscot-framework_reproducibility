Mapping cells through time, space, and beyond
=============================================
.. figure:: images/light_mode_logo.png
   :width: 400px
   :alt: The moscot toolbox.
   :align: center
   :figclass: center

   `moscot` is a toolbox for Optimal Transport (OT) applications in single-cell genomics. Moscot makes OT applications accessible to biologists by providing a unified framwork. Being built upon [OTT-JAX](https://ott-jax.readthedocs.io/en/latest/index.html), moscot leverages state-of-the-art advances in OT algorithms allowing for the application to large-scale single-cell datasets while enabling the incorporation of multiple modalities. While the documentation of the [software package](https://github.com/theislab/moscot) can be found on [moscot-tools.org](moscot-tools.org), this repository contains the data analyses perfomed for the _preprint.

Manuscript
----------
Our manuscript is available as a `preprint`_ on bioRxiv. 


Reproducibility
---------------
To ease reproducibility of our preprint results, the structure of this repository follows the one of the manuscript.
Each folder contains notebooks and scripts necessary to reproduce the corresponding analysis. 

Results
-------

.. csv-table::
   :header: "Application", "Folder path"

    moscot.time (Fig. 2), `notebooks/time/ <notebooks/time/>`__
    moscot.space (Fig. 3), `notebooks/space/ <notebooks/space/>`__
    moscot.spatiotemporal (Fig. 4), `notebooks/spatiotemporal/ <notebooks/spatiotemporal/>`__
    pancreas multiome (Fig. 5), `notebooks/pancreas/ <notebooks/spatiotemporal/>`__


.. _manuscript: TODO
