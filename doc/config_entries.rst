=================
The west.cfg File
=================

------------------
The [data] Section
------------------

The ``[data]`` section of ``west.cfg`` controls how data such as trajectory
connectivity and progress coordinate data is stored by WEST.

west_data_file (pathname)
  The name of the master WEST HDF5 file. (Default: ``west.h5``.)

west_data_file_driver (string)
  The name of the HDF5 low-level `file driver`_ used to access the WEST HDF5
  file. The default is to let h5py_ decide which driver to use. The only other
  potentially-useful value is ``family``, which will split the HDF5 file into
  2 GiB chunks.

  .. h5py: http://code.google.com/p/h5py/
  .. _`file driver`: http://h5py.alfven.org/docs-2.0/high/file.html#file-drivers

iter_prec (integer)
  The width of the iteration number field used in naming iteration groups in
  the HDF5 file. This defaults to 8, so that iteration groups are named like
  ``iter_00000001``. (This is done as a convenience so that iteration groups
  sort properly when viewing the WEST HDF5 file in hdfview_.)

  .. _hdfview: http://www.hdfgroup.org/hdf-java-html/hdfview/

aux_compression_threshold (integer)
  The minimum size of an auxiliary data set that will trigger automatic
  compression of the data set in the HDF5 file. (Default: 1 MiB.)

load_auxdata (boolean)
  By default, WEST does not load auxiliary data sets when loading segment
  data during the propagation/reweighting loop. If custom code requires
  auxiliary data to be present during propagation/reweighting, this parameter
  can be set to True to cause WEST to load auxiliary data into the ``data``
  dictionary of ``Segment`` objects. (Default: False.)
