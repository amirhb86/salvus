salvus - waveform propagator with benefits
==========================================

Hello and welcome to Salvus!

Impementing your own shit
^^^^^^^^^^^^^^^^^^^^^^^^^

Before running in and implementing your own stuff, its worthwhile to take a
quick tour of the important classes.  Depending on what you're trying to do,
you may only have to modify one, two, or all, of the different important
classes.

Problem
^^^^^^^

This is as close to a 'god' class as we get in Salvus. It has three main
functions: given a mesh, a reference element, a set of sources, and a model:

1.  Ensure that all of these different passed classes are compatible with one another.
2.  Set up the spatial discretization on each element, and apply the model parameters
3.  Run a time loop.




Contents:

.. toctree::
    :maxdepth: 2

    some_other_page
    api_index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
