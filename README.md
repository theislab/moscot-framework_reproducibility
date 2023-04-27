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


Code, tutorials and data
-------------------------
Under the hood,moslin is based on `moscot`_ to solve the optimal transport problem of mapping
lineage-traced cells across time points. Specifically, we implement moslin via the
`LineageClass`_ , we demonstrate a use case in our `tutorial`_ and we showcase
how to work with `tree distances`_ in an example. Downstream analysis, like
`visualizing`_ the inferred cell-cell transitions, is available via moscot's API.

Raw published data is available from the Gene Expression Omnibus (GEO) under accession codes:

- `c elegans`_: `GSE126954 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126954>`_.
- `zebrafish`_: `GSE159032  <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159032>`_.

Additionally, we simulated data using `LineageOT`_ and `TedSim`_. Processed data
is available on `figshare`_. To ease reproducibility, our data examples can
also be accessed through moscot's `dataset interface <https://moscot.readthedocs.io/en/latest/user.html#module-moscot.datasets>`_.

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
    moscot.spatiotemporal (Fig. 4), `notebooks/spatiotemporal/ <notebooks/spatiotemporal>`__
    pancreas multiome (Fig. 5), `notebooks/pancreas/ <notebooks/spatiotemporal>`__




.. _manuscript: TODO
.. _moscot: https://moscot-tools.org/
.. _preprint: https://www.biorxiv.org/content/10.1101/2023.04.14.536867v1
.. _LineageClass: https://moscot.readthedocs.io/en/latest/genapi/moscot.problems.time.LineageProblem.html
.. _tree distances: https://moscot.readthedocs.io/en/latest/notebooks/examples/problems/600_leaf_distance.html
.. _tutorial: https://moscot.readthedocs.io/en/latest/notebooks/tutorials/100_lineage.html
.. _LineageOT: https://doi.org/10.1038/s41467-021-25133-1
.. _TedSim: https://doi.org/10.1093/nar/gkac235
.. _c elegans: https://doi.org/10.1126/science.aax1971
.. _zebrafish: https://doi.org/10.1038/s41588-022-01129-5
.. _visualizing: https://moscot.readthedocs.io/en/latest/user.html#module-moscot.plotting
.. _figshare: https://doi.org/10.6084/m9.figshare.c.6533377.v1
