#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>  // For NumPy array support
#include <complex.h>            // For complex number support

// Declare Fortran subroutines (using the appropriate names if Fortran appends an underscore)
extern void permittivity_write_(int *im, double *mlist, int npoints_pml, int nx, int nz);
extern void permittivity_write_c_(int *im, double complex *mlist, int npoints_pml, int nx, int nz);
extern void attenuation_write_(int *im, double *mlist, int npoints_pml, int nx, int nz);
extern void stiffness_write_(int *im, double *mlist, int npoints_pml, int nx, int nz);
extern void seismic2_(int nx, int nz, double dx, double dz, int npoints_pml, int *src, int nstep, int single_output);
extern void seismic25_(int nx, int ny, int nz, double dx, double dy, double dz, int npoints_pml, int *src, int nstep, int single_output);
extern void electromag2_(int nx, int nz, double dx, double dz, int npoints_pml, int *src, int nstep, int single_output);
extern void electromag25_(int nx, int ny, int nz, double dx, double dy, double dz, int npoints_pml, int *src, int nstep, int single_output);
extern void electromag2c_(int nx, int nz, double dx, double dz, int npoints_pml, int *src, int nstep, int single_output);
extern void electromag25c_(int nx, int ny, int nz, double dx, double dy, double dz, int npoints_pml, int *src, int nstep, int single_output);

// =============================================================================
static PyObject *cpmlfdtd_permittivity_write(PyObject *self, PyObject *args) {
    PyArrayObject *im_obj, *mlist_obj;
    int npoints_pml, nx, nz;

    if (!PyArg_ParseTuple(args, "O!O!iii", 
                          &PyArray_Type, &im_obj, 
                          &PyArray_Type, &mlist_obj, 
                          &npoints_pml, &nx, &nz)) {
        return NULL;
    }

    if (PyArray_NDIM(im_obj) != 2 || PyArray_TYPE(im_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D integer array for 'im'");
        return NULL;
    }
    if (PyArray_NDIM(mlist_obj) != 2 || PyArray_TYPE(mlist_obj) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D double array for 'mlist'");
        return NULL;
    }

    int *im = (int *) PyArray_DATA(im_obj);
    double *mlist = (double *) PyArray_DATA(mlist_obj);

    permittivity_write_(im, mlist, npoints_pml, nx, nz);

    Py_RETURN_NONE;
}

// -----------------------------------------------------------------------------
// Wrapper for the Fortran subroutine permittivity_write_c
static PyObject *cpmlfdtd_permittivity_write_c(PyObject *self, PyObject *args) {
    PyArrayObject *im_obj, *mlist_obj;
    int npoints_pml, nx, nz;

    if (!PyArg_ParseTuple(args, "O!O!iii", 
                          &PyArray_Type, &im_obj, 
                          &PyArray_Type, &mlist_obj, 
                          &npoints_pml, &nx, &nz)) {
        return NULL;
    }

    if (PyArray_NDIM(im_obj) != 2 || PyArray_TYPE(im_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D integer array for 'im'");
        return NULL;
    }
    if (PyArray_NDIM(mlist_obj) != 2 || PyArray_TYPE(mlist_obj) != NPY_COMPLEX128) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D complex array for 'mlist'");
        return NULL;
    }

    int *im = (int *) PyArray_DATA(im_obj);
    double complex *mlist = (double complex *) PyArray_DATA(mlist_obj);

    permittivity_write_c_(im, mlist, npoints_pml, nx, nz);

    Py_RETURN_NONE;
}

// -----------------------------------------------------------------------------
// Wrapper for the Fortran subroutine attenuation_write
static PyObject *cpmlfdtd_attenuation_write(PyObject *self, PyObject *args) {
    PyArrayObject *im_obj, *mlist_obj;
    int npoints_pml, nx, nz;

    if (!PyArg_ParseTuple(args, "O!O!iii", 
                          &PyArray_Type, &im_obj, 
                          &PyArray_Type, &mlist_obj, 
                          &npoints_pml, &nx, &nz)) {
        return NULL;
    }

    if (PyArray_NDIM(im_obj) != 2 || PyArray_TYPE(im_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D integer array for 'im'");
        return NULL;
    }
    if (PyArray_NDIM(mlist_obj) != 2 || PyArray_TYPE(mlist_obj) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D double array for 'mlist'");
        return NULL;
    }

    int *im = (int *) PyArray_DATA(im_obj);
    double *mlist = (double *) PyArray_DATA(mlist_obj);

    attenuation_write_(im, mlist, npoints_pml, nx, nz);

    Py_RETURN_NONE;
}

// -----------------------------------------------------------------------------
// Wrapper for the Fortran subroutine stiffness_write
static PyObject *cpmlfdtd_stiffness_write(PyObject *self, PyObject *args) {
    PyArrayObject *im_obj, *mlist_obj;
    int npoints_pml, nx, nz;

    if (!PyArg_ParseTuple(args, "O!O!iii", 
                          &PyArray_Type, &im_obj, 
                          &PyArray_Type, &mlist_obj, 
                          &npoints_pml, &nx, &nz)) {
        return NULL;
    }

    if (PyArray_NDIM(im_obj) != 2 || PyArray_TYPE(im_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D integer array for 'im'");
        return NULL;
    }
    if (PyArray_NDIM(mlist_obj) != 2 || PyArray_TYPE(mlist_obj) != NPY_DOUBLE) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D double array for 'mlist'");
        return NULL;
    }

    int *im = (int *) PyArray_DATA(im_obj);
    double *mlist = (double *) PyArray_DATA(mlist_obj);

    stiffness_write_(im, mlist, npoints_pml, nx, nz);

    Py_RETURN_NONE;
}

// -----------------------------------------------------------------------------
// Wrapper for the Fortran subroutine seismic2
static PyObject *cpmlfdtd_seismic2(PyObject *self, PyObject *args) {
    int nx, nz, npoints_pml, nstep;
    double dx, dz;
    PyArrayObject *src_obj;
    int *src;
    int single_output = 1;

    // Parse the arguments from Python
    if (!PyArg_ParseTuple(args, "iidOi|p", 
                          &nx, &nz,      // Grid dimensions
                          &dx, &dz,      // Grid spacing
                          &src_obj,      // Source location array
                          &npoints_pml,  // PML thickness
                          &nstep,        // Number of time steps
                          &single_output)) {  // Optional: single precision output flag (default: True)
        return NULL;
    }

    // Ensure that src_obj is a 1D integer array
    if (PyArray_NDIM(src_obj) != 1 || PyArray_TYPE(src_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 1D integer array for 'src'");
        return NULL;
    }

    // Get the source array data
    src = (int *) PyArray_DATA(src_obj);

    // Call the Fortran subroutine
    seismic2_(nx, nz, dx, dz, npoints_pml, src, nstep, single_output);

    Py_RETURN_NONE;
}

// -----------------------------------------------------------------------------
// Wrapper for the Fortran subroutine seismic25
static PyObject *cpmlfdtd_seismic25(PyObject *self, PyObject *args) {
    int nx, ny, nz, npoints_pml, nstep;
    double dx, dy, dz;
    PyArrayObject *src_obj;
    int *src;
    int single_output = 1;  // Default to single precision output (True)

    // Parse the arguments from Python
    if (!PyArg_ParseTuple(args, "iiiidddOi|p",
                          &nx, &ny, &nz,    // Grid dimensions
                          &dx, &dy, &dz,    // Grid spacing
                          &src_obj,         // Source location array
                          &npoints_pml,     // PML thickness
                          &nstep,           // Number of time steps
                          &single_output)) {  // Optional: single precision output flag (default: True)
        return NULL;
    }

    // Ensure that src_obj is a 1D integer array
    if (PyArray_NDIM(src_obj) != 1 || PyArray_TYPE(src_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 1D integer array for 'src'");
        return NULL;
    }

    // Get the source array data
    src = (int *) PyArray_DATA(src_obj);

    // Call the Fortran subroutine
    seismic25_(nx, ny, nz, dx, dy, dz, npoints_pml, src, nstep, single_output);

    Py_RETURN_NONE;
}

// -----------------------------------------------------------------------------
// Wrapper for the Fortran subroutine electromag2
static PyObject *cpmlfdtd_electromag2(PyObject *self, PyObject *args) {
    int nx, nz, npoints_pml, nstep;
    double dx, dz;
    PyArrayObject *src_obj;
    int *src;
    int single_output = 1;

    // Parse the arguments from Python
    if (!PyArg_ParseTuple(args, "iidOi|p", 
                          &nx, &nz,      // Grid dimensions
                          &dx, &dz,      // Grid spacing
                          &src_obj,      // Source location array
                          &npoints_pml,  // PML thickness
                          &nstep,        // Number of time steps
                          &single_output)) {  // Optional: single precision output flag (default: True)
        return NULL;
    }

    // Ensure that src_obj is a 1D integer array
    if (PyArray_NDIM(src_obj) != 1 || PyArray_TYPE(src_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 1D integer array for 'src'");
        return NULL;
    }

    // Get the source array data
    src = (int *) PyArray_DATA(src_obj);

    // Call the Fortran subroutine
    electromag2_(nx, nz, dx, dz, npoints_pml, src, nstep, single_output);

    Py_RETURN_NONE;
}
// -----------------------------------------------------------------------------
// Wrapper for the Fortran subroutine electromag25
static PyObject *cpmlfdtd_electromag25(PyObject *self, PyObject *args) {
    int nx, ny, nz, npoints_pml, nstep;
    double dx, dy, dz;
    PyArrayObject *src_obj;
    int *src;
    int single_output = 1;  // Default to single precision output (True)

    // Parse the arguments from Python
    if (!PyArg_ParseTuple(args, "iiiidddOi|p",
                          &nx, &ny, &nz,    // Grid dimensions
                          &dx, &dy, &dz,    // Grid spacing
                          &src_obj,         // Source location array
                          &npoints_pml,     // PML thickness
                          &nstep,           // Number of time steps
                          &single_output)) {  // Optional: single precision output flag (default: True)
        return NULL;
    }

    // Ensure that src_obj is a 1D integer array
    if (PyArray_NDIM(src_obj) != 1 || PyArray_TYPE(src_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 1D integer array for 'src'");
        return NULL;
    }

    // Get the source array data
    src = (int *) PyArray_DATA(src_obj);

    // Call the Fortran subroutine
    seismic25_(nx, ny, nz, dx, dy, dz, npoints_pml, src, nstep, single_output);

    Py_RETURN_NONE;
}

// -----------------------------------------------------------------------------
// Wrapper for the Fortran subroutine electromag2c
static PyObject *cpmlfdtd_electromag2c(PyObject *self, PyObject *args) {
    int nx, nz, npoints_pml, nstep;
    double dx, dz;
    PyArrayObject *src_obj;
    int *src;
    int single_output = 1;  // Default to single precision output (True)

    // Parse the arguments from Python
    if (!PyArg_ParseTuple(args, "iidOi|p", 
                          &nx, &nz,      // Grid dimensions
                          &dx, &dz,      // Grid spacing
                          &src_obj,      // Source location array
                          &npoints_pml,  // PML thickness
                          &nstep,        // Number of time steps
                          &single_output)) {  // Optional: single precision output flag (default: True)
        return NULL;
    }

    // Ensure that src_obj is a 1D integer array
    if (PyArray_NDIM(src_obj) != 1 || PyArray_TYPE(src_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 1D integer array for 'src'");
        return NULL;
    }

    // Get the source array data
    src = (int *) PyArray_DATA(src_obj);

    // Call the Fortran subroutine
    electromag2c_(nx, nz, dx, dz, npoints_pml, src, nstep, single_output);

    Py_RETURN_NONE;
}

// -----------------------------------------------------------------------------
// Wrapper for the Fortran subroutine electromag25c
static PyObject *cpmlfdtd_electromag25c(PyObject *self, PyObject *args) {
    int nx, ny, nz, npoints_pml, nstep;
    double dx, dy, dz;
    PyArrayObject *src_obj;
    int *src;
    int single_output = 1;  // Default to single precision output (True)

    // Parse the arguments from Python
    if (!PyArg_ParseTuple(args, "iiiidddOi|p",
                          &nx, &ny, &nz,    // Grid dimensions
                          &dx, &dy, &dz,    // Grid spacing
                          &src_obj,         // Source location array
                          &npoints_pml,     // PML thickness
                          &nstep,           // Number of time steps
                          &single_output)) {  // Optional: single precision output flag (default: True)
        return NULL;
    }

    // Ensure that src_obj is a 1D integer array
    if (PyArray_NDIM(src_obj) != 1 || PyArray_TYPE(src_obj) != NPY_INT32) {
        PyErr_SetString(PyExc_TypeError, "Expected a 1D integer array for 'src'");
        return NULL;
    }

    // Get the source array data
    src = (int *) PyArray_DATA(src_obj);

    // Call the Fortran subroutine
    electromag25c_(nx, ny, nz, dx, dy, dz, npoints_pml, src, nstep, single_output);

    Py_RETURN_NONE;
}
// =============================================================================
// Define the methods for the Python module
static PyMethodDef CpmlfdtdMethods[] = {
    {"permittivity_write", cpmlfdtd_permittivity_write, METH_VARARGS, "Write permittivity"},
    {"permittivity_write_c", cpmlfdtd_permittivity_write_c, METH_VARARGS, "Write permittivity with complex values"},
    {"attenuation_write", cpmlfdtd_attenuation_write, METH_VARARGS, "Write attenuation"},
    {"seismic2", cpmlfdtd_seismic2, METH_VARARGS, "Run 2D seismic simulation"},
    {"electromag2", cpmlfdtd_electromag2, METH_VARARGS, "Run 2D electromagnetic simulation"},
    {"seismic25", cpmlfdtd_seismic25, METH_VARARGS, "Run 2.5D seismic simulation"},
    {"electromag25", cpmlfdtd_electromag25, METH_VARARGS, "Run 2.5D electromagnetic simulation"},
    {"electromag2c", cpmlfdtd_electromag2c, METH_VARARGS, "Run 2D electromagnetic simulation with complex values"},
    {"electromag25c", cpmlfdtd_electromag25c, METH_VARARGS, "Run 2.5D electromagnetic simulation with complex values"},
    {"stiffness_write", cpmlfdtd_stiffness_write, METH_VARARGS, "Write stiffness"},
    {NULL, NULL, 0, NULL}  // Sentinel
};


// Define the Python module
static struct PyModuleDef cpmlfdtdmodule = {
    PyModuleDef_HEAD_INIT,
    "cpmlfdtd",  // Name of the module
    NULL,        // Module documentation (can be NULL)
    -1,          // Size of per-interpreter state of the module, or -1 if the module keeps state in global variables
    CpmlfdtdMethods
};

// Module initialization function
PyMODINIT_FUNC PyInit_cpmlfdtd(void) {
    import_array();  // Initialize NumPy
    return PyModule_Create(&cpmlfdtdmodule);
}