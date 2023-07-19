// Not required for MSVC since the code is included below
#if defined(_WIN32) && !defined(_MSC_VER)
double lgamma(double x);
#endif

// lgamma not included in ANSI standard and so not available in MSVC
#if defined(_MSC_VER)
double lgamma(double z) {
    double c[7];
    double x,y ,tmp, ser, v;
    int i;

    if (z<=0) return 0;

    c[0]=2.5066282746310005;
    c[1]=76.18009172947146;
    c[2]=-86.50532032941677;
    c[3]=24.01409824083091;
    c[4]=-1.231739572450155;
    c[5]=0.1208650973866179e-2;
    c[6]=-0.5395239384953e-5;

    x   = z;
    y   = x;
    tmp = x+5.5;
    tmp = (x+0.5)*log(tmp)-tmp;
    ser = 1.000000000190015;
    for (i=1; i<7; i++) {
        y   += 1.0;
        ser += c[i]/y;
        }
    v = tmp+log(c[0]*ser/x);
    return v;
    }
#endif

static PyObject *fact_wrap(PyObject *self,PyObject *args){
  int ok = 0, n=0;
  ok = PyArg_ParseTuple(args,"i",&n);
  if (!ok) return NULL;
  return Py_BuildValue("i",fact(n));
}
static PyObject *fact2_wrap(PyObject *self,PyObject *args){
  int ok = 0, n=0;
  ok = PyArg_ParseTuple(args,"i",&n);
  if (!ok) return NULL;
  return Py_BuildValue("i",fact2(n));
}
static PyObject *dist2_wrap(PyObject *self,PyObject *args){
  int ok = 0;
  double x1,y1,z1,x2,y2,z2;
  PyObject *A, *B;
  ok = PyArg_ParseTuple(args,"OO",&A,&B);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(A,"ddd",&x1,&y1,&z1);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(B,"ddd",&x2,&y2,&z2);
  if (!ok) return NULL;
  return Py_BuildValue("d",dist2(x1,y1,z1,x2,y2,z2));
}
static PyObject *dist_wrap(PyObject *self,PyObject *args){
  int ok = 0;
  double x1,y1,z1,x2,y2,z2;
  PyObject *A, *B;
  ok = PyArg_ParseTuple(args,"OO",&A,&B);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(A,"ddd",&x1,&y1,&z1);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(B,"ddd",&x2,&y2,&z2);
  if (!ok) return NULL;
  return Py_BuildValue("d",dist(x1,y1,z1,x2,y2,z2));
}
static PyObject *binomial_wrap(PyObject *self,PyObject *args){
  int ok = 0, ia=0, ib=0;
  ok = PyArg_ParseTuple(args,"ii",&ia,&ib);
  if (!ok) return NULL;
  return Py_BuildValue("i",binomial(ia,ib));
}
static PyObject *binomial_prefactor_wrap(PyObject *self,PyObject *args){
  int ok = 0, s=0, ia=0, ib=0;
  double xpa=0.,xpb=0.;
  ok = PyArg_ParseTuple(args,"iiidd",&s,&ia,&ib,&xpa,&xpb);
  if (!ok) return NULL;
  return Py_BuildValue("i",binomial_prefactor(s,ia,ib,xpa,xpb));
}
static PyObject *Fgamma_wrap(PyObject *self,PyObject *args){
  int ok = 0;
  double m=0.,x=0.;
  ok = PyArg_ParseTuple(args,"dd",&m,&x);
  if (!ok) return NULL;
  return Py_BuildValue("d",Fgamma(m,x));
}
static PyObject *ijkl2intindex_wrap(PyObject *self,PyObject *args){
  int ok = 0,i,j,k,l;
  ok = PyArg_ParseTuple(args,"iiii",&i,&j,&k,&l);
  if (!ok) return NULL;
  return Py_BuildValue("i",ijkl2intindex(i,j,k,l));
}
static PyObject *fB_wrap(PyObject *self,PyObject *args){
  int ok = 0,i,l1,l2,r;
  double px,ax,bx,g;
  ok = PyArg_ParseTuple(args,"iiidddid",&i,&l1,&l2,&px,&ax,&bx,&r,&g);
  if (!ok) return NULL;
  return Py_BuildValue("d",fB(i,l1,l2,px,ax,bx,r,g));
}
static PyObject *fact_ratio2_wrap(PyObject *self,PyObject *args){
  int ok = 0,a,b;
  ok = PyArg_ParseTuple(args,"ii",&a,&b);
  if (!ok) return NULL;
  return Py_BuildValue("i",fact_ratio2(a,b));
}

static PyObject *contr_coulomb_wrap(PyObject *self,PyObject *args){
  int ok=0;
  double xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd;
  int lena,lenb,lenc,lend;
  PyObject *aexps_obj,*acoefs_obj,*anorms_obj,
    *bexps_obj,*bcoefs_obj,*bnorms_obj,
    *cexps_obj,*ccoefs_obj,*cnorms_obj,
    *dexps_obj,*dcoefs_obj,*dnorms_obj,
    *xyza_obj,*lmna_obj,*xyzb_obj,*lmnb_obj,
    *xyzc_obj,*lmnc_obj,*xyzd_obj,*lmnd_obj;
  double *aexps,*acoefs,*anorms,
    *bexps,*bcoefs,*bnorms,
    *cexps,*ccoefs,*cnorms,
    *dexps,*dcoefs,*dnorms;
  int i;
  double Jij=0; /* return value */

  ok = PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOO",
			&aexps_obj,&acoefs_obj,&anorms_obj,&xyza_obj,&lmna_obj,
			&bexps_obj,&bcoefs_obj,&bnorms_obj,&xyzb_obj,&lmnb_obj,
			&cexps_obj,&ccoefs_obj,&cnorms_obj,&xyzc_obj,&lmnc_obj,
			&dexps_obj,&dcoefs_obj,&dnorms_obj,&xyzd_obj,&lmnd_obj);
  if (!ok) return NULL;


  ok=PyArg_ParseTuple(xyza_obj,"ddd",&xa,&ya,&za);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(xyzb_obj,"ddd",&xb,&yb,&zb);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(xyzc_obj,"ddd",&xc,&yc,&zc);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(xyzd_obj,"ddd",&xd,&yd,&zd);
  if (!ok) return NULL;

  ok=PyArg_ParseTuple(lmna_obj,"iii",&la,&ma,&na);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(lmnb_obj,"iii",&lb,&mb,&nb);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(lmnc_obj,"iii",&lc,&mc,&nc);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(lmnd_obj,"iii",&ld,&md,&nd);
  if (!ok) return NULL;

  /* Test that each is a sequence: */
  if (!PySequence_Check(aexps_obj)) return NULL;
  if (!PySequence_Check(acoefs_obj)) return NULL;
  if (!PySequence_Check(anorms_obj)) return NULL;
  if (!PySequence_Check(bexps_obj)) return NULL;
  if (!PySequence_Check(bcoefs_obj)) return NULL;
  if (!PySequence_Check(bnorms_obj)) return NULL;
  if (!PySequence_Check(cexps_obj)) return NULL;
  if (!PySequence_Check(ccoefs_obj)) return NULL;
  if (!PySequence_Check(cnorms_obj)) return NULL;
  if (!PySequence_Check(dexps_obj)) return NULL;
  if (!PySequence_Check(dcoefs_obj)) return NULL;
  if (!PySequence_Check(dnorms_obj)) return NULL;

  /* Get the length of each sequence */
  lena = PySequence_Size(aexps_obj);
  if (lena<0) return NULL;
  if (lena != PySequence_Size(acoefs_obj)) return NULL;
  if (lena != PySequence_Size(anorms_obj)) return NULL;
  lenb = PySequence_Size(bexps_obj);
  if (lenb<0) return NULL;
  if (lenb != PySequence_Size(bcoefs_obj)) return NULL;
  if (lenb != PySequence_Size(bnorms_obj)) return NULL;
  lenc = PySequence_Size(cexps_obj);
  if (lenc<0) return NULL;
  if (lenc != PySequence_Size(ccoefs_obj)) return NULL;
  if (lenc != PySequence_Size(cnorms_obj)) return NULL;
  lend = PySequence_Size(dexps_obj);
  if (lend<0) return NULL;
  if (lend != PySequence_Size(dcoefs_obj)) return NULL;
  if (lend != PySequence_Size(dnorms_obj)) return NULL;

  /* Allocate the space for each array */
  if (lena+lenb+lenc+lend > 4*MAX_PRIMS_PER_CONT) return NULL;
  aexps = work;
  acoefs = aexps + lena;
  anorms = acoefs + lena;
  bexps = anorms + lena;
  bcoefs = bexps + lenb;
  bnorms = bcoefs + lenb;
  cexps = bnorms + lenb;
  ccoefs = cexps + lenc;
  cnorms = ccoefs + lenc;
  dexps = cnorms + lenc;
  dcoefs = dexps + lend;
  dnorms = dcoefs + lend;

  /* Unpack all of the lengths: */
  for (i=0; i<lena; i++){
    aexps[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(aexps_obj,i));
    acoefs[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(acoefs_obj,i));
    anorms[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(anorms_obj,i));
  }
  for (i=0; i<lenb; i++){
    bexps[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(bexps_obj,i));
    bcoefs[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(bcoefs_obj,i));
    bnorms[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(bnorms_obj,i));
  }
  for (i=0; i<lenc; i++){
    cexps[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(cexps_obj,i));
    ccoefs[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(ccoefs_obj,i));
    cnorms[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(cnorms_obj,i));
  }
  for (i=0; i<lend; i++){
    dexps[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(dexps_obj,i));
    dcoefs[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(dcoefs_obj,i));
    dnorms[i] = PyFloat_AS_DOUBLE(PySequence_GetItem(dnorms_obj,i));
  }
  Jij = contr_coulomb(lena,aexps,acoefs,anorms,xa,ya,za,la,ma,na,
		      lenb,bexps,bcoefs,bnorms,xb,yb,zb,lb,mb,nb,
		      lenc,cexps,ccoefs,cnorms,xc,yc,zc,lc,mc,nc,
		      lend,dexps,dcoefs,dnorms,xd,yd,zd,ld,md,nd);

  return Py_BuildValue("d", Jij);
}

static PyObject *contr_nuke_vec_wrap(PyObject *self,PyObject *args){

  /* This turned out to be slower than multiple calls to
     nuclear_attraction_vec. I'm leaving the code in place
     for reference.
  */
  PyObject *aexps, *acoefs, *anorms, *aorigin, *apowers,
    *bexps, *bcoefs, *bnorms, *borigin, *bpowers,
    *xc, *yc, *zc, *wc, *qc;
  int ok;
  double xa,ya,za,xb,yb,zb;
  int la,ma,na,lb,mb,nb;
  int nprima, nprimb, ncenters;
  int i,j,k;
  double anormi,aexpi,acoefi, bnormj,bexpj,bcoefj,xck,yck,zck,wck,qck;
  double incr=0, Vnij=0;


  ok = PyArg_ParseTuple(args,"OOOOOOOOOOOOOOO",
			&aexps,&acoefs,&anorms,&aorigin,&apowers,
			&bexps,&bcoefs,&bnorms,&borigin,&bpowers,
			&xc, &yc, &zc, &wc, &qc);
  if (!ok) return NULL;

  ok = PyArg_ParseTuple(aorigin,"ddd",&xa,&ya,&za);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(borigin,"ddd",&xb,&yb,&zb);
  if (!ok) return NULL;

  ok = PyArg_ParseTuple(apowers,"iii",&la,&ma,&na);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(bpowers,"iii",&lb,&mb,&nb);
  if (!ok) return NULL;

  nprima = PySequence_Size(aexps);
  if (nprima<0) return NULL;
  if (nprima != PySequence_Size(acoefs)) return NULL;
  if (nprima != PySequence_Size(anorms)) return NULL;

  nprimb = PySequence_Size(bexps);
  if (nprimb<0) return NULL;
  if (nprimb != PySequence_Size(bcoefs)) return NULL;
  if (nprimb != PySequence_Size(bnorms)) return NULL;

  ncenters = PySequence_Size(xc);
  if (ncenters < 0) return NULL;
  if (ncenters != PySequence_Size(yc)) return NULL;
  if (ncenters != PySequence_Size(zc)) return NULL;
  if (ncenters != PySequence_Size(wc)) return NULL;
  if (ncenters != PySequence_Size(qc)) return NULL;

  for (k=0; k<ncenters; k++){
    xck = PyFloat_AS_DOUBLE(PySequence_GetItem(xc,k));
    yck = PyFloat_AS_DOUBLE(PySequence_GetItem(yc,k));
    zck = PyFloat_AS_DOUBLE(PySequence_GetItem(zc,k));
    wck = PyFloat_AS_DOUBLE(PySequence_GetItem(wc,k));
    qck = PyFloat_AS_DOUBLE(PySequence_GetItem(qc,k));
    for (i=0; i<nprima; i++){
      anormi = PyFloat_AS_DOUBLE(PySequence_GetItem(anorms,i));
      aexpi = PyFloat_AS_DOUBLE(PySequence_GetItem(aexps,i));
      acoefi = PyFloat_AS_DOUBLE(PySequence_GetItem(acoefs,i));
      for (j=0; j<nprimb; j++){
	bnormj = PyFloat_AS_DOUBLE(PySequence_GetItem(bnorms,j));
	bexpj = PyFloat_AS_DOUBLE(PySequence_GetItem(bexps,j));
	bcoefj = PyFloat_AS_DOUBLE(PySequence_GetItem(bcoefs,j));
	incr = nuclear_attraction(xa,ya,za,anormi,la,ma,na,aexpi,
				  xb,yb,zb,bnormj,lb,mb,nb,bexpj,
				  xck,yck,zck);
	Vnij += acoefi*bcoefj*wck*qck*incr;
      }
    }
  }
  return Py_BuildValue("d",Vnij);

}

static PyObject *coulomb_repulsion_wrap(PyObject *self,PyObject *args){
  int ok=0;
  double norma,alphaa,normb,alphab,normc,alphac,normd,alphad,
    xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd;
  int la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd;
  PyObject *A,*B,*C,*D,*powa,*powb,*powc,*powd;

  ok=PyArg_ParseTuple(args,"OdOdOdOdOdOdOdOd",&A,&norma,&powa,&alphaa,
		      &B,&normb,&powb,&alphab,&C,&normc,&powc,&alphac,
		      &D,&normd,&powd,&alphad);
  if (!ok) return NULL;

  ok=PyArg_ParseTuple(A,"ddd",&xa,&ya,&za);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(B,"ddd",&xb,&yb,&zb);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(C,"ddd",&xc,&yc,&zc);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(D,"ddd",&xd,&yd,&zd);
  if (!ok) return NULL;

  ok=PyArg_ParseTuple(powa,"iii",&la,&ma,&na);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(powb,"iii",&lb,&mb,&nb);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(powc,"iii",&lc,&mc,&nc);
  if (!ok) return NULL;
  ok=PyArg_ParseTuple(powd,"iii",&ld,&md,&nd);
  if (!ok) return NULL;

  return Py_BuildValue("d",
    coulomb_repulsion(xa,ya,za,norma,la,ma,na,alphaa,
		      xb,yb,zb,normb,lb,mb,nb,alphab,
		      xc,yc,zc,normc,lc,mc,nc,alphac,
		      xd,yd,zd,normd,ld,md,nd,alphad));
}

static PyObject *kinetic_wrap(PyObject *self,PyObject *args){
  int ok=0,l1,m1,n1,l2,m2,n2;
  double xa,ya,za,xb,yb,zb,alpha1,alpha2;
  PyObject *A,*B,*powa,*powb;

  ok = PyArg_ParseTuple(args,"dOOdOO",&alpha1,&powa,&A,&alpha2,
			&powb,&B);

  if (!ok) return NULL;

  ok = PyArg_ParseTuple(powa,"iii",&l1,&m1,&n1);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(powb,"iii",&l2,&m2,&n2);
  if (!ok) return NULL;

  ok = PyArg_ParseTuple(A,"ddd",&xa,&ya,&za);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(B,"ddd",&xb,&yb,&zb);
  if (!ok) return NULL;

  return Py_BuildValue("d",
		       kinetic(alpha1,l1,m1,n1,xa,ya,za,
			       alpha2,l2,m2,n2,xb,yb,zb));
}

static PyObject *overlap_wrap(PyObject *self,PyObject *args){
  int ok=0,l1,m1,n1,l2,m2,n2;
  double xa,ya,za,xb,yb,zb,alpha1,alpha2;
  PyObject *A,*B,*powa,*powb;

  ok = PyArg_ParseTuple(args,"dOOdOO",&alpha1,&powa,&A,&alpha2,
			&powb,&B);
  if (!ok) return NULL;

  ok = PyArg_ParseTuple(powa,"iii",&l1,&m1,&n1);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(powb,"iii",&l2,&m2,&n2);
  if (!ok) return NULL;

  ok = PyArg_ParseTuple(A,"ddd",&xa,&ya,&za);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(B,"ddd",&xb,&yb,&zb);
  if (!ok) return NULL;

  return Py_BuildValue("d",
		       overlap(alpha1,l1,m1,n1,xa,ya,za,
				      alpha2,l2,m2,n2,xb,yb,zb));
}

static PyObject *nuclear_attraction_wrap(PyObject *self,PyObject *args){
  int ok=0,l1,m1,n1,l2,m2,n2;
  double x1,y1,z1,x2,y2,z2,x3,y3,z3,norm1,alpha1,norm2,alpha2;
  PyObject *A,*B,*C,*powa,*powb;

  ok = PyArg_ParseTuple(args,"OdOdOdOdO",&A,&norm1,&powa,&alpha1,
			&B,&norm2,&powb,&alpha2,&C);
  if (!ok) return NULL;

  ok = PyArg_ParseTuple(A,"ddd",&x1,&y1,&z1);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(B,"ddd",&x2,&y2,&z2);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(C,"ddd",&x3,&y3,&z3);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(powa,"iii",&l1,&m1,&n1);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(powb,"iii",&l2,&m2,&n2);
  if (!ok) return NULL;

  return Py_BuildValue("d",
     nuclear_attraction(x1,y1,z1,norm1,l1,m1,n1,alpha1,
			x2,y2,z2,norm2,l2,m2,n2,alpha2,
			x3,y3,z3));
}

static PyObject *three_center_1D_wrap(PyObject *self,PyObject *args){

  double xi,xj,xk,alphai,alphaj,alphak;
  int ai,aj,ak;
  int ok=0;

  ok = PyArg_ParseTuple(args,"diddiddid",&xi,&ai,&alphai,
			&xj,&aj,&alphaj,&xk,&ak,&alphak);
  if (!ok) return NULL;
  return Py_BuildValue("d",
		       three_center_1D(xi,ai,alphai,
				       xj,aj,alphaj,xk,ak,alphak));
}

static PyObject *nuclear_attraction_vec_wrap(PyObject *self,PyObject *args){
  int ok=0,l1,m1,n1,l2,m2,n2;
  double x1,y1,z1,x2,y2,z2,norm1,alpha1,norm2,alpha2;
  PyObject *A,*B,*xc_obj, *yc_obj, *zc_obj, *wc_obj, *qc_obj, 
    *powa, *powb;
  double retval = 0,wc,qc,xc,yc,zc;
  int veclength = 0, i;

  ok = PyArg_ParseTuple(args,"OdOdOdOdOOOOO",&A,&norm1,&powa,&alpha1,
			&B,&norm2,&powb,&alpha2,
			&xc_obj, &yc_obj, &zc_obj, &wc_obj, &qc_obj);
  if (!ok) return NULL;

  ok = PyArg_ParseTuple(A,"ddd",&x1,&y1,&z1);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(B,"ddd",&x2,&y2,&z2);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(powa,"iii",&l1,&m1,&n1);
  if (!ok) return NULL;
  ok = PyArg_ParseTuple(powb,"iii",&l2,&m2,&n2);
  if (!ok) return NULL;

  if (!PySequence_Check(xc_obj)) return NULL;
  if (!PySequence_Check(yc_obj)) return NULL;
  if (!PySequence_Check(zc_obj)) return NULL;
  if (!PySequence_Check(wc_obj)) return NULL;
  if (!PySequence_Check(qc_obj)) return NULL;


  veclength = PySequence_Size(xc_obj);
  if (veclength<0) return NULL;
  if (veclength != PySequence_Size(yc_obj)) return NULL;
  if (veclength != PySequence_Size(zc_obj)) return NULL;
  if (veclength != PySequence_Size(wc_obj)) return NULL;
  if (veclength != PySequence_Size(qc_obj)) return NULL;


  for (i=0; i<veclength; i++){
    xc = PyFloat_AS_DOUBLE(PySequence_GetItem(xc_obj,i));
    yc = PyFloat_AS_DOUBLE(PySequence_GetItem(yc_obj,i));
    zc = PyFloat_AS_DOUBLE(PySequence_GetItem(zc_obj,i));
    wc = PyFloat_AS_DOUBLE(PySequence_GetItem(wc_obj,i));
    qc = PyFloat_AS_DOUBLE(PySequence_GetItem(qc_obj,i));

    retval += wc*qc*nuclear_attraction(x1,y1,z1,norm1,l1,m1,n1,alpha1,
				     x2,y2,z2,norm2,l2,m2,n2,alpha2,
				     xc,yc,zc);
  }

  return Py_BuildValue("d",retval);
}

/* Python interface */
static PyMethodDef cints_methods[] = {
  {"fact",fact_wrap,METH_VARARGS},
  {"fact2",fact2_wrap,METH_VARARGS},
  {"dist2",dist2_wrap,METH_VARARGS},
  {"dist",dist_wrap,METH_VARARGS},
  {"binomial",binomial_wrap,METH_VARARGS},
  {"binomial_prefactor",binomial_prefactor_wrap,METH_VARARGS},
  {"Fgamma",Fgamma_wrap,METH_VARARGS},
  {"ijkl2intindex",ijkl2intindex_wrap,METH_VARARGS},
  {"fB",fB_wrap,METH_VARARGS},
  {"fact_ratio2",fact_ratio2_wrap,METH_VARARGS},
  {"contr_coulomb",contr_coulomb_wrap,METH_VARARGS},
  {"coulomb_repulsion",coulomb_repulsion_wrap,METH_VARARGS},
  {"kinetic",kinetic_wrap,METH_VARARGS},
  {"overlap",overlap_wrap,METH_VARARGS},
  {"nuclear_attraction",nuclear_attraction_wrap,METH_VARARGS},
  {"nuclear_attraction_vec",nuclear_attraction_vec_wrap,METH_VARARGS},
  {"contr_nuke_vec",contr_nuke_vec_wrap,METH_VARARGS},
  {"three_center_1D",three_center_1D_wrap,METH_VARARGS},
  {NULL,NULL} /* Sentinel */
};

static void module_init(char* name)
{
  (void) Py_InitModule(name,cints_methods);
}

#if defined(_WIN32)
__declspec(dllexport)
#endif
#if defined(PYQUANTE_FULLY_QUALIFIED_MODULE_NAME)
void initpyquante_cints_ext(){module_init("pyquante_cints_ext");}
#else
void initcints(){module_init("cints");}
#endif
