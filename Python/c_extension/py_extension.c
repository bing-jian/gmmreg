/**
 * $Author: bing.jian $
 * $Date: 2008-11-06 15:45:11 -0500 (Thu, 06 Nov 2008) $
 * $Revision: 98 $
 */

#include <Python.h>

/*
Reference:
  Writing a C extension to NumPy
  http://numpy.scipy.org/numpydoc/numpy-13.html
Note the location of arrayobject.h has been changed! See
http://projects.scipy.org/pipermail/numpy-tickets/2006-October/000344.html
In [13]: numpy.version.version
Out[13]: '1.1.1'
In [14]: numpy.get_include()
Out[14]: 'C:\\DevTools\\Python25\\lib\\site-packages\\numpy\\core\\include'
*/
#include "numpy/arrayobject.h"


#define IND2(a,i,j) \
    *((double *)(a->data + i*a->strides[0] + j*a->strides[1]))

#include "cvgmi_API.h"

static PyObject *
py_squared_distance_matrix(PyObject *self, PyObject *args)
{
    int m,n,d;
    double *dist;
    PyObject *A, *B, *g;
    PyArrayObject *arrayA, *arrayB, *arrayg;
    PyArrayObject *array_dist;
    /*
    double scale, result;
    PyObject *list;
    */
    int out_dim[2];
    /* send Python arrays to C */
    if (!PyArg_ParseTuple(args, "OOO",  &A, &B, &g))
    {
        return NULL;
    }

    /* printf("getting input success. \n"); */
    arrayA = (PyArrayObject *) PyArray_ContiguousFromObject(A, PyArray_DOUBLE, 1, 2);
    arrayB = (PyArrayObject *) PyArray_ContiguousFromObject(B, PyArray_DOUBLE, 1, 2);
    arrayg = (PyArrayObject *) PyArray_ContiguousFromObject(g, PyArray_DOUBLE, 1, 2);
    /* printf("converting success!\n"); */

    if (arrayA->nd > 2 || arrayA->descr->type_num != PyArray_DOUBLE) {
        PyErr_SetString(PyExc_ValueError,
        "array must be two-dimensional and of type float");
        return NULL;
    }
    /*printf("checking input success!\n");*/

    m = (arrayA->dimensions)[0];
    n = (arrayB->dimensions)[0];
    if (arrayA->nd>1)
        d = (arrayA->dimensions)[1];
    else
        d = 1;
    /* printf("m=%d,n=%d,d=%d\n",m,n,d); */
    dist = (double *) malloc(m*n*sizeof(double));

    /* call function */
    squared_distance_matrix((double*)(arrayA->data), (double*)(arrayB->data), (double*)(arrayg->data), m, n, d, dist);
    out_dim[0] = m;
    out_dim[1] = n;
    
    /* PyArray_FromDimsAndData() deprecated, use PyArray_SimpleNewFromData()
    http://blog.enthought.com/?p=62
    array_dist = (PyArrayObject*) PyArray_FromDimsAndData(2,out_dim,PyArray_DOUBLE, (char*)dist);
    */
    array_dist =  PyArray_SimpleNewFromData(2,out_dim,NPY_DOUBLE,dist);

    if (array_dist == NULL){
        printf("creating %dx%d array failed\n", out_dim[0],out_dim[1]);
        return NULL;
    }
    /*free(dist);*/

    /* send the result back to Python */
    Py_DECREF(arrayA);
    Py_DECREF(arrayB);
    Py_DECREF(arrayg);

    return PyArray_Return(array_dist);
    //return PyFloat_FromDouble(result);

    /*return Py_BuildValue("d", result);*/
}

static PyObject *
py_gauss_transform(PyObject *self, PyObject *args)
{
    int m,n,dim;
    double scale, result;
    double *grad;
    PyObject *A, *B;
    PyArrayObject *arrayA, *arrayB;
    PyArrayObject *arrayGrad;
    PyObject *list;

    /* send Python arrays to C */
    if (!PyArg_ParseTuple(args, "OOd",  &A, &B, &scale))
    {
        return NULL;
    }

    /* printf("getting input success. \n"); */
    arrayA = (PyArrayObject *) PyArray_ContiguousFromObject(A, PyArray_DOUBLE, 1, 2);
    arrayB = (PyArrayObject *) PyArray_ContiguousFromObject(B, PyArray_DOUBLE, 1, 2);

    /* printf("converting success!\n"); */

    if (arrayA->nd > 2 || arrayA->descr->type_num != PyArray_DOUBLE) {
        PyErr_SetString(PyExc_ValueError,
        "array must be two-dimensional and of type float");
        return NULL;
    }
    /* printf("checking input success!\n"); */

    m = (arrayA->dimensions)[0];
    n = (arrayB->dimensions)[0];
    if (arrayA->nd>1)
        dim = (arrayA->dimensions)[1];
    else
        dim = 1;
    /* printf("m=%d,n=%d,dim=%d\n",m,n,dim); */
    grad = (double *) malloc(m*dim*sizeof(double));

    /* call function */
    result = GaussTransform((double*)(arrayA->data), (double*)(arrayB->data), m, n, dim, scale, grad);

    /* PyArray_FromDimsAndData() deprecated, use PyArray_SimpleNewFromData()
    http://blog.enthought.com/?p=62
    arrayGrad = (PyArrayObject*) PyArray_FromDimsAndData(2,arrayA->dimensions,PyArray_DOUBLE, (char*)grad);
    */
    arrayGrad = PyArray_SimpleNewFromData(2,arrayA->dimensions,NPY_DOUBLE,grad);
    if (arrayGrad == NULL){
        printf("creating %dx%d array failed\n", arrayA->dimensions[0],arrayA->dimensions[1]);
        return NULL;
    }
    /* free(grad); */

    /* send the result back to Python */
    Py_DECREF(arrayA);
    Py_DECREF(arrayB);


    //Build a list; send it back to interpreter
    list = PyList_New(0);
    // Check the API documentation for meaning of
    // return values.
    if(PyList_Append(list, PyFloat_FromDouble(result)) != 0)
    {
     // set exception context, raise (return 0)
        return 0;
    }
    if(PyList_Append(list, PyArray_Return(arrayGrad)) != 0)
    {
     // set exception context, raise (return 0)
        return 0;
    }

    return list;
    //return PyArray_Return(arrayGrad);
    //return PyFloat_FromDouble(result);

    /*return Py_BuildValue("d", result);*/
}



static PyMethodDef pyMethods[] = {
    {"squared_distance_matrix",  py_squared_distance_matrix, METH_VARARGS,
     "Compute the squared distance matrix."},
    {"gauss_transform",  py_gauss_transform, METH_VARARGS,
     "Compute the Gauss Transform."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
init_extension(void)
{
    (void) Py_InitModule("_extension", pyMethods);
    import_array();
}

