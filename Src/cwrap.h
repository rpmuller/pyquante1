/* My routines */
#ifdef _MSC_VER
double lgamma(double);
#endif

/* Wrappers */
static PyObject *fact_wrap(PyObject *self,PyObject *args);
static PyObject *fact2_wrap(PyObject *self,PyObject *args);
static PyObject *dist2_wrap(PyObject *self,PyObject *args);
static PyObject *dist_wrap(PyObject *self,PyObject *args);
static PyObject *dist_wrap(PyObject *self,PyObject *args);
static PyObject *binomial_prefactor_wrap(PyObject *self,PyObject *args);
static PyObject *Fgamma_wrap(PyObject *self,PyObject *args);
static PyObject *ijkl2intindex_wrap(PyObject *self,PyObject *args);
static PyObject *fB_wrap(PyObject *self,PyObject *args);
static PyObject *fact_ratio2_wrap(PyObject *self,PyObject *args);
static PyObject *contr_coulomb_wrap(PyObject *self,PyObject *args);
static PyObject *coulomb_repulsion_wrap(PyObject *self,PyObject *args);
static PyObject *kinetic_wrap(PyObject *self,PyObject *args);
static PyObject *overlap_wrap(PyObject *self,PyObject *args);
static PyObject *nuclear_attraction_wrap(PyObject *self,PyObject *args);
static PyObject *nuclear_attraction_vec_wrap(PyObject *self,PyObject *args);
static PyObject *three_center_1D_wrap(PyObject *self,PyObject *args);

