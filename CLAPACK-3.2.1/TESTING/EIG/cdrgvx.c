/* cdrgvx.f -- translated by f2c (version 20061008).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"
#include "blaswrap.h"

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static complex c_b11 = {1.f,0.f};
static integer c__5 = 5;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static integer c__3 = 3;
static integer c__6 = 6;
static integer c__4 = 4;

/* Subroutine */ int cdrgvx_(integer *nsize, real *thresh, integer *nin, 
	integer *nout, complex *a, integer *lda, complex *b, complex *ai, 
	complex *bi, complex *alpha, complex *beta, complex *vl, complex *vr, 
	integer *ilo, integer *ihi, real *lscale, real *rscale, real *s, real 
	*stru, real *dif, real *diftru, complex *work, integer *lwork, real *
	rwork, integer *iwork, integer *liwork, real *result, logical *bwork, 
	integer *info)
{
    /* Format strings */
    static char fmt_9999[] = "(\002 CDRGVX: \002,a,\002 returned INFO=\002,i"
	    "6,\002.\002,/9x,\002N=\002,i6,\002, JTYPE=\002,i6,\002)\002)";
    static char fmt_9998[] = "(\002 CDRGVX: \002,a,\002 Eigenvectors from"
	    " \002,a,\002 incorrectly \002,\002normalized.\002,/\002 Bits of "
	    "error=\002,0p,g10.3,\002,\002,9x,\002N=\002,i6,\002, JTYPE=\002,"
	    "i6,\002, IWA=\002,i5,\002, IWB=\002,i5,\002, IWX=\002,i5,\002, I"
	    "WY=\002,i5)";
    static char fmt_9997[] = "(/1x,a3,\002 -- Complex Expert Eigenvalue/vect"
	    "or\002,\002 problem driver\002)";
    static char fmt_9995[] = "(\002 Matrix types: \002,/)";
    static char fmt_9994[] = "(\002 TYPE 1: Da is diagonal, Db is identity,"
	    " \002,/\002     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1) \002,/"
	    "\002     YH and X are left and right eigenvectors. \002,/)";
    static char fmt_9993[] = "(\002 TYPE 2: Da is quasi-diagonal, Db is iden"
	    "tity, \002,/\002     A = Y^(-H) Da X^(-1), B = Y^(-H) Db X^(-1)"
	    " \002,/\002     YH and X are left and right eigenvectors. \002,/)"
	    ;
    static char fmt_9992[] = "(/\002 Tests performed:  \002,/4x,\002 a is al"
	    "pha, b is beta, l is a left eigenvector, \002,/4x,\002 r is a ri"
	    "ght eigenvector and \002,a,\002 means \002,a,\002.\002,/\002 1 ="
	    " max | ( b A - a B )\002,a,\002 l | / const.\002,/\002 2 = max |"
	    " ( b A - a B ) r | / const.\002,/\002 3 = max ( Sest/Stru, Stru/"
	    "Sest ) \002,\002 over all eigenvalues\002,/\002 4 = max( DIFest/"
	    "DIFtru, DIFtru/DIFest ) \002,\002 over the 1st and 5th eigenvect"
	    "ors\002,/)";
    static char fmt_9991[] = "(\002 Type=\002,i2,\002,\002,\002 IWA=\002,i2"
	    ",\002, IWB=\002,i2,\002, IWX=\002,i2,\002, IWY=\002,i2,\002, res"
	    "ult \002,i2,\002 is\002,0p,f8.2)";
    static char fmt_9990[] = "(\002 Type=\002,i2,\002,\002,\002 IWA=\002,i2"
	    ",\002, IWB=\002,i2,\002, IWX=\002,i2,\002, IWY=\002,i2,\002, res"
	    "ult \002,i2,\002 is\002,1p,e10.3)";
    static char fmt_9987[] = "(\002 CDRGVX: \002,a,\002 returned INFO=\002,i"
	    "6,\002.\002,/9x,\002N=\002,i6,\002, Input example #\002,i2,\002"
	    ")\002)";
    static char fmt_9986[] = "(\002 CDRGVX: \002,a,\002 Eigenvectors from"
	    " \002,a,\002 incorrectly \002,\002normalized.\002,/\002 Bits of "
	    "error=\002,0p,g10.3,\002,\002,9x,\002N=\002,i6,\002, Input Examp"
	    "le #\002,i2,\002)\002)";
    static char fmt_9996[] = "(\002Input Example\002)";
    static char fmt_9989[] = "(\002 Input example #\002,i2,\002, matrix orde"
	    "r=\002,i4,\002,\002,\002 result \002,i2,\002 is\002,0p,f8.2)";
    static char fmt_9988[] = "(\002 Input example #\002,i2,\002, matrix orde"
	    "r=\002,i4,\002,\002,\002 result \002,i2,\002 is\002,1p,e10.3)";

    /* System generated locals */
    integer a_dim1, a_offset, ai_dim1, ai_offset, b_dim1, b_offset, bi_dim1, 
	    bi_offset, vl_dim1, vl_offset, vr_dim1, vr_offset, i__1, i__2;
    real r__1, r__2, r__3, r__4;
    complex q__1;

    /* Builtin functions */
    double sqrt(doublereal);
    void c_div(complex *, complex *, complex *);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void),
	     s_rsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsle(void);

    /* Local variables */
    integer i__, j, n, iwa, iwb;
    real ulp;
    integer iwx, iwy, nmax;
    extern /* Subroutine */ int cget52_(logical *, integer *, complex *, 
	    integer *, complex *, integer *, complex *, integer *, complex *, 
	    complex *, complex *, real *, real *);
    integer linfo;
    real anorm, bnorm;
    integer nerrs;
    extern /* Subroutine */ int clatm6_(integer *, integer *, complex *, 
	    integer *, complex *, complex *, integer *, complex *, integer *, 
	    complex *, complex *, complex *, complex *, real *, real *);
    real ratio1, ratio2, thrsh2;
    extern doublereal clange_(char *, integer *, integer *, complex *, 
	    integer *, real *), slamch_(char *);
    extern /* Subroutine */ int clacpy_(char *, integer *, integer *, complex 
	    *, integer *, complex *, integer *), xerbla_(char *, 
	    integer *);
    real abnorm;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *);
    extern /* Subroutine */ int alasvm_(char *, integer *, integer *, integer 
	    *, integer *), cggevx_(char *, char *, char *, char *, 
	    integer *, complex *, integer *, complex *, integer *, complex *, 
	    complex *, complex *, integer *, complex *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, 
	    complex *, integer *, real *, integer *, logical *, integer *);
    complex weight[5];
    integer minwrk, maxwrk, iptype;
    real ulpinv;
    integer nptknt, ntestt;

    /* Fortran I/O blocks */
    static cilist io___20 = { 0, 0, 0, fmt_9999, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_9998, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_9998, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_9997, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_9995, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_9994, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_9993, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_9992, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_9991, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_9990, 0 };
    static cilist io___35 = { 0, 0, 1, 0, 0 };
    static cilist io___36 = { 0, 0, 0, 0, 0 };
    static cilist io___37 = { 0, 0, 0, 0, 0 };
    static cilist io___38 = { 0, 0, 0, 0, 0 };
    static cilist io___39 = { 0, 0, 0, 0, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_9987, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_9986, 0 };
    static cilist io___42 = { 0, 0, 0, fmt_9986, 0 };
    static cilist io___43 = { 0, 0, 0, fmt_9997, 0 };
    static cilist io___44 = { 0, 0, 0, fmt_9996, 0 };
    static cilist io___45 = { 0, 0, 0, fmt_9992, 0 };
    static cilist io___46 = { 0, 0, 0, fmt_9989, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_9988, 0 };



/*  -- LAPACK test routine (version 3.1) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
/*     November 2006 */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CDRGVX checks the nonsymmetric generalized eigenvalue problem */
/*  expert driver CGGEVX. */

/*  CGGEVX computes the generalized eigenvalues, (optionally) the left */
/*  and/or right eigenvectors, (optionally) computes a balancing */
/*  transformation to improve the conditioning, and (optionally) */
/*  reciprocal condition numbers for the eigenvalues and eigenvectors. */

/*  When CDRGVX is called with NSIZE > 0, two types of test matrix pairs */
/*  are generated by the subroutine SLATM6 and test the driver CGGEVX. */
/*  The test matrices have the known exact condition numbers for */
/*  eigenvalues. For the condition numbers of the eigenvectors */
/*  corresponding the first and last eigenvalues are also know */
/*  ``exactly'' (see CLATM6). */
/*  For each matrix pair, the following tests will be performed and */
/*  compared with the threshhold THRESH. */

/*  (1) max over all left eigenvalue/-vector pairs (beta/alpha,l) of */

/*     | l**H * (beta A - alpha B) | / ( ulp max( |beta A|, |alpha B| ) ) */

/*      where l**H is the conjugate tranpose of l. */

/*  (2) max over all right eigenvalue/-vector pairs (beta/alpha,r) of */

/*        | (beta A - alpha B) r | / ( ulp max( |beta A|, |alpha B| ) ) */

/*  (3) The condition number S(i) of eigenvalues computed by CGGEVX */
/*      differs less than a factor THRESH from the exact S(i) (see */
/*      CLATM6). */

/*  (4) DIF(i) computed by CTGSNA differs less than a factor 10*THRESH */
/*      from the exact value (for the 1st and 5th vectors only). */

/*  Test Matrices */
/*  ============= */

/*  Two kinds of test matrix pairs */
/*           (A, B) = inverse(YH) * (Da, Db) * inverse(X) */
/*  are used in the tests: */

/*  1: Da = 1+a   0    0    0    0    Db = 1   0   0   0   0 */
/*           0   2+a   0    0    0         0   1   0   0   0 */
/*           0    0   3+a   0    0         0   0   1   0   0 */
/*           0    0    0   4+a   0         0   0   0   1   0 */
/*           0    0    0    0   5+a ,      0   0   0   0   1 , and */

/*  2: Da =  1   -1    0    0    0    Db = 1   0   0   0   0 */
/*           1    1    0    0    0         0   1   0   0   0 */
/*           0    0    1    0    0         0   0   1   0   0 */
/*           0    0    0   1+a  1+b        0   0   0   1   0 */
/*           0    0    0  -1-b  1+a ,      0   0   0   0   1 . */

/*  In both cases the same inverse(YH) and inverse(X) are used to compute */
/*  (A, B), giving the exact eigenvectors to (A,B) as (YH, X): */

/*  YH:  =  1    0   -y    y   -y    X =  1   0  -x  -x   x */
/*          0    1   -y    y   -y         0   1   x  -x  -x */
/*          0    0    1    0    0         0   0   1   0   0 */
/*          0    0    0    1    0         0   0   0   1   0 */
/*          0    0    0    0    1,        0   0   0   0   1 , where */

/*  a, b, x and y will have all values independently of each other from */
/*  { sqrt(sqrt(ULP)),  0.1,  1,  10,  1/sqrt(sqrt(ULP)) }. */

/*  Arguments */
/*  ========= */

/*  NSIZE   (input) INTEGER */
/*          The number of sizes of matrices to use.  NSIZE must be at */
/*          least zero. If it is zero, no randomly generated matrices */
/*          are tested, but any test matrices read from NIN will be */
/*          tested.  If it is not zero, then N = 5. */

/*  THRESH  (input) REAL */
/*          A test will count as "failed" if the "error", computed as */
/*          described above, exceeds THRESH.  Note that the error */
/*          is scaled to be O(1), so THRESH should be a reasonably */
/*          small multiple of 1, e.g., 10 or 100.  In particular, */
/*          it should not depend on the precision (single vs. double) */
/*          or the size of the matrix.  It must be at least zero. */

/*  NIN     (input) INTEGER */
/*          The FORTRAN unit number for reading in the data file of */
/*          problems to solve. */

/*  NOUT    (input) INTEGER */
/*          The FORTRAN unit number for printing out error messages */
/*          (e.g., if a routine returns IINFO not equal to 0.) */

/*  A       (workspace) COMPLEX array, dimension (LDA, NSIZE) */
/*          Used to hold the matrix whose eigenvalues are to be */
/*          computed.  On exit, A contains the last matrix actually used. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of A, B, AI, BI, Ao, and Bo. */
/*          It must be at least 1 and at least NSIZE. */

/*  B       (workspace) COMPLEX array, dimension (LDA, NSIZE) */
/*          Used to hold the matrix whose eigenvalues are to be */
/*          computed.  On exit, B contains the last matrix actually used. */

/*  AI      (workspace) COMPLEX array, dimension (LDA, NSIZE) */
/*          Copy of A, modified by CGGEVX. */

/*  BI      (workspace) COMPLEX array, dimension (LDA, NSIZE) */
/*          Copy of B, modified by CGGEVX. */

/*  ALPHA   (workspace) COMPLEX array, dimension (NSIZE) */
/*  BETA    (workspace) COMPLEX array, dimension (NSIZE) */
/*          On exit, ALPHA/BETA are the eigenvalues. */

/*  VL      (workspace) COMPLEX array, dimension (LDA, NSIZE) */
/*          VL holds the left eigenvectors computed by CGGEVX. */

/*  VR      (workspace) COMPLEX array, dimension (LDA, NSIZE) */
/*          VR holds the right eigenvectors computed by CGGEVX. */

/*  ILO     (output/workspace) INTEGER */

/*  IHI     (output/workspace) INTEGER */

/*  LSCALE  (output/workspace) REAL array, dimension (N) */

/*  RSCALE  (output/workspace) REAL array, dimension (N) */

/*  S       (output/workspace) REAL array, dimension (N) */

/*  STRU    (output/workspace) REAL array, dimension (N) */

/*  DIF     (output/workspace) REAL array, dimension (N) */

/*  DIFTRU  (output/workspace) REAL array, dimension (N) */

/*  WORK    (workspace) COMPLEX array, dimension (LWORK) */

/*  LWORK   (input) INTEGER */
/*          Leading dimension of WORK.  LWORK >= 2*N*N + 2*N */

/*  RWORK   (workspace) REAL array, dimension (6*N) */

/*  IWORK   (workspace) INTEGER array, dimension (LIWORK) */

/*  LIWORK  (input) INTEGER */
/*          Leading dimension of IWORK.  LIWORK >= N+2. */

/*  RESULT  (output/workspace) REAL array, dimension (4) */

/*  BWORK   (workspace) LOGICAL array, dimension (N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */
/*          > 0:  A routine returned an error code. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check for errors */

    /* Parameter adjustments */
    vr_dim1 = *lda;
    vr_offset = 1 + vr_dim1;
    vr -= vr_offset;
    vl_dim1 = *lda;
    vl_offset = 1 + vl_dim1;
    vl -= vl_offset;
    bi_dim1 = *lda;
    bi_offset = 1 + bi_dim1;
    bi -= bi_offset;
    ai_dim1 = *lda;
    ai_offset = 1 + ai_dim1;
    ai -= ai_offset;
    b_dim1 = *lda;
    b_offset = 1 + b_dim1;
    b -= b_offset;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --alpha;
    --beta;
    --lscale;
    --rscale;
    --s;
    --stru;
    --dif;
    --diftru;
    --work;
    --rwork;
    --iwork;
    --result;
    --bwork;

    /* Function Body */
    *info = 0;

    nmax = 5;

    if (*nsize < 0) {
	*info = -1;
    } else if (*thresh < 0.f) {
	*info = -2;
    } else if (*nin <= 0) {
	*info = -3;
    } else if (*nout <= 0) {
	*info = -4;
    } else if (*lda < 1 || *lda < nmax) {
	*info = -6;
    } else if (*liwork < nmax + 2) {
	*info = -26;
    }

/*     Compute workspace */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

    minwrk = 1;
    if (*info == 0 && *lwork >= 1) {
	minwrk = (nmax << 1) * (nmax + 1);
	maxwrk = nmax * (ilaenv_(&c__1, "CGEQRF", " ", &nmax, &c__1, &nmax, &
		c__0) + 1);
/* Computing MAX */
	i__1 = maxwrk, i__2 = (nmax << 1) * (nmax + 1);
	maxwrk = max(i__1,i__2);
	work[1].r = (real) maxwrk, work[1].i = 0.f;
    }

    if (*lwork < minwrk) {
	*info = -23;
    }

    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("CDRGVX", &i__1);
	return 0;
    }

    n = 5;
    ulp = slamch_("P");
    ulpinv = 1.f / ulp;
    thrsh2 = *thresh * 10.f;
    nerrs = 0;
    nptknt = 0;
    ntestt = 0;

    if (*nsize == 0) {
	goto L90;
    }

/*     Parameters used for generating test matrices. */

    r__1 = sqrt(sqrt(ulp));
    q__1.r = r__1, q__1.i = 0.f;
    weight[0].r = q__1.r, weight[0].i = q__1.i;
    weight[1].r = .1f, weight[1].i = 0.f;
    weight[2].r = 1.f, weight[2].i = 0.f;
    c_div(&q__1, &c_b11, &weight[1]);
    weight[3].r = q__1.r, weight[3].i = q__1.i;
    c_div(&q__1, &c_b11, weight);
    weight[4].r = q__1.r, weight[4].i = q__1.i;

    for (iptype = 1; iptype <= 2; ++iptype) {
	for (iwa = 1; iwa <= 5; ++iwa) {
	    for (iwb = 1; iwb <= 5; ++iwb) {
		for (iwx = 1; iwx <= 5; ++iwx) {
		    for (iwy = 1; iwy <= 5; ++iwy) {

/*                    generated a pair of test matrix */

			clatm6_(&iptype, &c__5, &a[a_offset], lda, &b[
				b_offset], &vr[vr_offset], lda, &vl[vl_offset]
, lda, &weight[iwa - 1], &weight[iwb - 1], &
				weight[iwx - 1], &weight[iwy - 1], &stru[1], &
				diftru[1]);

/*                    Compute eigenvalues/eigenvectors of (A, B). */
/*                    Compute eigenvalue/eigenvector condition numbers */
/*                    using computed eigenvectors. */

			clacpy_("F", &n, &n, &a[a_offset], lda, &ai[ai_offset]
, lda);
			clacpy_("F", &n, &n, &b[b_offset], lda, &bi[bi_offset]
, lda);

			cggevx_("N", "V", "V", "B", &n, &ai[ai_offset], lda, &
				bi[bi_offset], lda, &alpha[1], &beta[1], &vl[
				vl_offset], lda, &vr[vr_offset], lda, ilo, 
				ihi, &lscale[1], &rscale[1], &anorm, &bnorm, &
				s[1], &dif[1], &work[1], lwork, &rwork[1], &
				iwork[1], &bwork[1], &linfo);
			if (linfo != 0) {
			    io___20.ciunit = *nout;
			    s_wsfe(&io___20);
			    do_fio(&c__1, "CGGEVX", (ftnlen)6);
			    do_fio(&c__1, (char *)&linfo, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer))
				    ;
			    do_fio(&c__1, (char *)&iptype, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwa, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwb, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwx, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwy, (ftnlen)sizeof(
				    integer));
			    e_wsfe();
			    goto L30;
			}

/*                    Compute the norm(A, B) */

			clacpy_("Full", &n, &n, &ai[ai_offset], lda, &work[1], 
				 &n);
			clacpy_("Full", &n, &n, &bi[bi_offset], lda, &work[n *
				 n + 1], &n);
			i__1 = n << 1;
			abnorm = clange_("Fro", &n, &i__1, &work[1], &n, &
				rwork[1]);

/*                    Tests (1) and (2) */

			result[1] = 0.f;
			cget52_(&c_true, &n, &a[a_offset], lda, &b[b_offset], 
				lda, &vl[vl_offset], lda, &alpha[1], &beta[1], 
				 &work[1], &rwork[1], &result[1]);
			if (result[2] > *thresh) {
			    io___22.ciunit = *nout;
			    s_wsfe(&io___22);
			    do_fio(&c__1, "Left", (ftnlen)4);
			    do_fio(&c__1, "CGGEVX", (ftnlen)6);
			    do_fio(&c__1, (char *)&result[2], (ftnlen)sizeof(
				    real));
			    do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer))
				    ;
			    do_fio(&c__1, (char *)&iptype, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwa, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwb, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwx, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwy, (ftnlen)sizeof(
				    integer));
			    e_wsfe();
			}

			result[2] = 0.f;
			cget52_(&c_false, &n, &a[a_offset], lda, &b[b_offset], 
				 lda, &vr[vr_offset], lda, &alpha[1], &beta[1]
, &work[1], &rwork[1], &result[2]);
			if (result[3] > *thresh) {
			    io___23.ciunit = *nout;
			    s_wsfe(&io___23);
			    do_fio(&c__1, "Right", (ftnlen)5);
			    do_fio(&c__1, "CGGEVX", (ftnlen)6);
			    do_fio(&c__1, (char *)&result[3], (ftnlen)sizeof(
				    real));
			    do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer))
				    ;
			    do_fio(&c__1, (char *)&iptype, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwa, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwb, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwx, (ftnlen)sizeof(
				    integer));
			    do_fio(&c__1, (char *)&iwy, (ftnlen)sizeof(
				    integer));
			    e_wsfe();
			}

/*                    Test (3) */

			result[3] = 0.f;
			i__1 = n;
			for (i__ = 1; i__ <= i__1; ++i__) {
			    if (s[i__] == 0.f) {
				if (stru[i__] > abnorm * ulp) {
				    result[3] = ulpinv;
				}
			    } else if (stru[i__] == 0.f) {
				if (s[i__] > abnorm * ulp) {
				    result[3] = ulpinv;
				}
			    } else {
/* Computing MAX */
				r__3 = (r__1 = stru[i__] / s[i__], dabs(r__1))
					, r__4 = (r__2 = s[i__] / stru[i__], 
					dabs(r__2));
				rwork[i__] = dmax(r__3,r__4);
/* Computing MAX */
				r__1 = result[3], r__2 = rwork[i__];
				result[3] = dmax(r__1,r__2);
			    }
/* L10: */
			}

/*                    Test (4) */

			result[4] = 0.f;
			if (dif[1] == 0.f) {
			    if (diftru[1] > abnorm * ulp) {
				result[4] = ulpinv;
			    }
			} else if (diftru[1] == 0.f) {
			    if (dif[1] > abnorm * ulp) {
				result[4] = ulpinv;
			    }
			} else if (dif[5] == 0.f) {
			    if (diftru[5] > abnorm * ulp) {
				result[4] = ulpinv;
			    }
			} else if (diftru[5] == 0.f) {
			    if (dif[5] > abnorm * ulp) {
				result[4] = ulpinv;
			    }
			} else {
/* Computing MAX */
			    r__3 = (r__1 = diftru[1] / dif[1], dabs(r__1)), 
				    r__4 = (r__2 = dif[1] / diftru[1], dabs(
				    r__2));
			    ratio1 = dmax(r__3,r__4);
/* Computing MAX */
			    r__3 = (r__1 = diftru[5] / dif[5], dabs(r__1)), 
				    r__4 = (r__2 = dif[5] / diftru[5], dabs(
				    r__2));
			    ratio2 = dmax(r__3,r__4);
			    result[4] = dmax(ratio1,ratio2);
			}

			ntestt += 4;

/*                    Print out tests which fail. */

			for (j = 1; j <= 4; ++j) {
			    if (result[j] >= thrsh2 && j >= 4 || result[j] >= 
				    *thresh && j <= 3) {

/*                       If this is the first test to fail, */
/*                       print a header to the data file. */

				if (nerrs == 0) {
				    io___28.ciunit = *nout;
				    s_wsfe(&io___28);
				    do_fio(&c__1, "CXV", (ftnlen)3);
				    e_wsfe();

/*                          Print out messages for built-in examples */

/*                          Matrix types */

				    io___29.ciunit = *nout;
				    s_wsfe(&io___29);
				    e_wsfe();
				    io___30.ciunit = *nout;
				    s_wsfe(&io___30);
				    e_wsfe();
				    io___31.ciunit = *nout;
				    s_wsfe(&io___31);
				    e_wsfe();

/*                          Tests performed */

				    io___32.ciunit = *nout;
				    s_wsfe(&io___32);
				    do_fio(&c__1, "'", (ftnlen)1);
				    do_fio(&c__1, "transpose", (ftnlen)9);
				    do_fio(&c__1, "'", (ftnlen)1);
				    e_wsfe();

				}
				++nerrs;
				if (result[j] < 1e4f) {
				    io___33.ciunit = *nout;
				    s_wsfe(&io___33);
				    do_fio(&c__1, (char *)&iptype, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&iwa, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&iwb, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&iwx, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&iwy, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(
					    integer));
				    do_fio(&c__1, (char *)&result[j], (ftnlen)
					    sizeof(real));
				    e_wsfe();
				} else {
				    io___34.ciunit = *nout;
				    s_wsfe(&io___34);
				    do_fio(&c__1, (char *)&iptype, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&iwa, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&iwb, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&iwx, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&iwy, (ftnlen)
					    sizeof(integer));
				    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(
					    integer));
				    do_fio(&c__1, (char *)&result[j], (ftnlen)
					    sizeof(real));
				    e_wsfe();
				}
			    }
/* L20: */
			}

L30:

/* L40: */
			;
		    }
/* L50: */
		}
/* L60: */
	    }
/* L70: */
	}
/* L80: */
    }

    goto L150;

L90:

/*     Read in data from file to check accuracy of condition estimation */
/*     Read input data until N=0 */

    io___35.ciunit = *nin;
    i__1 = s_rsle(&io___35);
    if (i__1 != 0) {
	goto L150;
    }
    i__1 = do_lio(&c__3, &c__1, (char *)&n, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L150;
    }
    i__1 = e_rsle();
    if (i__1 != 0) {
	goto L150;
    }
    if (n == 0) {
	goto L150;
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___36.ciunit = *nin;
	s_rsle(&io___36);
	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    do_lio(&c__6, &c__1, (char *)&a[i__ + j * a_dim1], (ftnlen)sizeof(
		    complex));
	}
	e_rsle();
/* L100: */
    }
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___37.ciunit = *nin;
	s_rsle(&io___37);
	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    do_lio(&c__6, &c__1, (char *)&b[i__ + j * b_dim1], (ftnlen)sizeof(
		    complex));
	}
	e_rsle();
/* L110: */
    }
    io___38.ciunit = *nin;
    s_rsle(&io___38);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__4, &c__1, (char *)&stru[i__], (ftnlen)sizeof(real));
    }
    e_rsle();
    io___39.ciunit = *nin;
    s_rsle(&io___39);
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__4, &c__1, (char *)&diftru[i__], (ftnlen)sizeof(real));
    }
    e_rsle();

    ++nptknt;

/*     Compute eigenvalues/eigenvectors of (A, B). */
/*     Compute eigenvalue/eigenvector condition numbers */
/*     using computed eigenvectors. */

    clacpy_("F", &n, &n, &a[a_offset], lda, &ai[ai_offset], lda);
    clacpy_("F", &n, &n, &b[b_offset], lda, &bi[bi_offset], lda);

    cggevx_("N", "V", "V", "B", &n, &ai[ai_offset], lda, &bi[bi_offset], lda, 
	    &alpha[1], &beta[1], &vl[vl_offset], lda, &vr[vr_offset], lda, 
	    ilo, ihi, &lscale[1], &rscale[1], &anorm, &bnorm, &s[1], &dif[1], 
	    &work[1], lwork, &rwork[1], &iwork[1], &bwork[1], &linfo);

    if (linfo != 0) {
	io___40.ciunit = *nout;
	s_wsfe(&io___40);
	do_fio(&c__1, "CGGEVX", (ftnlen)6);
	do_fio(&c__1, (char *)&linfo, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nptknt, (ftnlen)sizeof(integer));
	e_wsfe();
	goto L140;
    }

/*     Compute the norm(A, B) */

    clacpy_("Full", &n, &n, &ai[ai_offset], lda, &work[1], &n);
    clacpy_("Full", &n, &n, &bi[bi_offset], lda, &work[n * n + 1], &n);
    i__1 = n << 1;
    abnorm = clange_("Fro", &n, &i__1, &work[1], &n, &rwork[1]);

/*     Tests (1) and (2) */

    result[1] = 0.f;
    cget52_(&c_true, &n, &a[a_offset], lda, &b[b_offset], lda, &vl[vl_offset], 
	     lda, &alpha[1], &beta[1], &work[1], &rwork[1], &result[1]);
    if (result[2] > *thresh) {
	io___41.ciunit = *nout;
	s_wsfe(&io___41);
	do_fio(&c__1, "Left", (ftnlen)4);
	do_fio(&c__1, "CGGEVX", (ftnlen)6);
	do_fio(&c__1, (char *)&result[2], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nptknt, (ftnlen)sizeof(integer));
	e_wsfe();
    }

    result[2] = 0.f;
    cget52_(&c_false, &n, &a[a_offset], lda, &b[b_offset], lda, &vr[vr_offset]
, lda, &alpha[1], &beta[1], &work[1], &rwork[1], &result[2]);
    if (result[3] > *thresh) {
	io___42.ciunit = *nout;
	s_wsfe(&io___42);
	do_fio(&c__1, "Right", (ftnlen)5);
	do_fio(&c__1, "CGGEVX", (ftnlen)6);
	do_fio(&c__1, (char *)&result[3], (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nptknt, (ftnlen)sizeof(integer));
	e_wsfe();
    }

/*     Test (3) */

    result[3] = 0.f;
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s[i__] == 0.f) {
	    if (stru[i__] > abnorm * ulp) {
		result[3] = ulpinv;
	    }
	} else if (stru[i__] == 0.f) {
	    if (s[i__] > abnorm * ulp) {
		result[3] = ulpinv;
	    }
	} else {
/* Computing MAX */
	    r__3 = (r__1 = stru[i__] / s[i__], dabs(r__1)), r__4 = (r__2 = s[
		    i__] / stru[i__], dabs(r__2));
	    rwork[i__] = dmax(r__3,r__4);
/* Computing MAX */
	    r__1 = result[3], r__2 = rwork[i__];
	    result[3] = dmax(r__1,r__2);
	}
/* L120: */
    }

/*     Test (4) */

    result[4] = 0.f;
    if (dif[1] == 0.f) {
	if (diftru[1] > abnorm * ulp) {
	    result[4] = ulpinv;
	}
    } else if (diftru[1] == 0.f) {
	if (dif[1] > abnorm * ulp) {
	    result[4] = ulpinv;
	}
    } else if (dif[5] == 0.f) {
	if (diftru[5] > abnorm * ulp) {
	    result[4] = ulpinv;
	}
    } else if (diftru[5] == 0.f) {
	if (dif[5] > abnorm * ulp) {
	    result[4] = ulpinv;
	}
    } else {
/* Computing MAX */
	r__3 = (r__1 = diftru[1] / dif[1], dabs(r__1)), r__4 = (r__2 = dif[1] 
		/ diftru[1], dabs(r__2));
	ratio1 = dmax(r__3,r__4);
/* Computing MAX */
	r__3 = (r__1 = diftru[5] / dif[5], dabs(r__1)), r__4 = (r__2 = dif[5] 
		/ diftru[5], dabs(r__2));
	ratio2 = dmax(r__3,r__4);
	result[4] = dmax(ratio1,ratio2);
    }

    ntestt += 4;

/*     Print out tests which fail. */

    for (j = 1; j <= 4; ++j) {
	if (result[j] >= thrsh2) {

/*           If this is the first test to fail, */
/*           print a header to the data file. */

	    if (nerrs == 0) {
		io___43.ciunit = *nout;
		s_wsfe(&io___43);
		do_fio(&c__1, "CXV", (ftnlen)3);
		e_wsfe();

/*              Print out messages for built-in examples */

/*              Matrix types */

		io___44.ciunit = *nout;
		s_wsfe(&io___44);
		e_wsfe();

/*              Tests performed */

		io___45.ciunit = *nout;
		s_wsfe(&io___45);
		do_fio(&c__1, "'", (ftnlen)1);
		do_fio(&c__1, "transpose", (ftnlen)9);
		do_fio(&c__1, "'", (ftnlen)1);
		e_wsfe();

	    }
	    ++nerrs;
	    if (result[j] < 1e4f) {
		io___46.ciunit = *nout;
		s_wsfe(&io___46);
		do_fio(&c__1, (char *)&nptknt, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&result[j], (ftnlen)sizeof(real));
		e_wsfe();
	    } else {
		io___47.ciunit = *nout;
		s_wsfe(&io___47);
		do_fio(&c__1, (char *)&nptknt, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&result[j], (ftnlen)sizeof(real));
		e_wsfe();
	    }
	}
/* L130: */
    }

L140:

    goto L90;
L150:

/*     Summary */

    alasvm_("CXV", nout, &nerrs, &ntestt, &c__0);

    work[1].r = (real) maxwrk, work[1].i = 0.f;

    return 0;















/*     End of CDRGVX */

} /* cdrgvx_ */
