Installation
============

Installation of iGFA can be done using pip:

.. code-block:: console

   (.venv) $ pip install igfa

In addition, iGFA relies on IPOPT solver for the NLP optimization.

.. note::
   IPOPT's default installation uses the MUMPS linear solver (which iGFA uses). To configure more powerful solvers, refer to `cyipopt documentation <https://cyipopt.readthedocs.io/en/stable/install.html>`_


Installing IPOPT for Linux
--------------------------

For Linux, IPOPT can be installed using conda:

.. code-block:: console

   (.venv) $ conda install -c conda-forge cyipopt

Installing IPOPT for Mac/Windows
--------------------------------

There are issues with installing the latest version of IPOPT on Windows/Mac. However, users can make use of the 3.11.1 version installable with conda as:
.. code-block:: console

   (.venv) $ conda install -c conda-forge ipopt==3.11.1

Testing Installation of IPOPT
-----------------------------
The following command should list the version of IPOPT that is installed. If this doesn't work, the IPOPT istallation has been unsuccessful.

.. code-block:: console

   (.venv) $ ipopt -v


