Analysing STORM data: ``storm`` accessors
=========================================

Importing anything from `xstorm` registers these accessors in your Python
environment. The accessor methods are designed to be called on Datasets opened
using :py:func:`open_stormdataset` which performs initialisation and creates
necessary metadata. The methods can be accessed for a Dataset ``ds`` through
``ds.storm``, e.g.

.. code-block:: python

    from xstorm import open_stormdataset
    ds = open_stormdataset("boutdata.nc")
    mfp = ds.storm.electron_mfp()

.. autoclass:: xstorm.StormDatasetAccessor
    :members:

.. autoclass:: xstorm.StormDataArrayAccessor
    :members:
