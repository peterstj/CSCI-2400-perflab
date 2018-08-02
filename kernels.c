/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"

/*
 * Please fill in the following student struct:
 */
const team_t team = {
    "tepe5782@colorado.edu", //Replace this with your email address.
    ""                   //Replace this with your partner's email address. Leave blank if working alone.
};

/***************
 * FLIP KERNEL
 ***************/

/******************************************************
 * Your different versions of the flip kernel go here
 ******************************************************/
char version_one_descr[] = "version one";
void version_one(int dim, pixel *src, pixel *dst)
{
    int i, j;
    for (i = 0; i < dim; i++){
        for (j = 0; j < dim; j++){
            dst[RIDX(i, j, dim)].red = src[RIDX(j, i, dim)].red;
	    dst[RIDX(i, j, dim)].green = src[RIDX(j, i, dim)].green;
	    dst[RIDX(i, j, dim)].blue = src[RIDX(j, i, dim)].blue;
        }
    }
    for (i = 0; i < dim; i++){
	for (j = 0; j < dim; j++){
	    dst[RIDX(i, j, dim)].red = src[RIDX(i, dim-1-j, dim)].red;
	    dst[RIDX(i, j, dim)].green = src[RIDX(i, dim-1-j, dim)].green;
	    dst[RIDX(i, j, dim)].blue = src[RIDX(i, dim-1-j, dim)].blue;
	}
    }
}

char v2_descr[] = "version two";
void v2(int dim, pixel *src, pixel *dst)
{
    int i, j, a;
    for (i = 0; i < dim; i++){
	for (j=0; j < dim; j++){
	    a = dim-1-i;
	    dst[RIDX(j, a, dim)].red = src[RIDX(i, j, dim)].red;
	    dst[RIDX(j, a, dim)].green = src[RIDX(i, j, dim)].green;
	    dst[RIDX(j, a, dim)].blue = src[RIDX(i, j, dim)].blue;
	}
    }
}

char v3_descr[] = "version three";
void v3(int dim, pixel *src, pixel *dst)
{
    int i, j, a;
    for (i = 0 ; i < dim; i++){
	for (j=0; j < dim; j++){
	    a = dim-1-i;
	    dst[RIDX(i, j, dim)].red = src[RIDX(j, a, dim)].red;
	    dst[RIDX(i, j, dim)].green = src[RIDX(j, a, dim)].green;
	    dst[RIDX(i, j, dim)].blue = src[RIDX(j, a, dim)].blue;
	}
    }
}
char v4_descr[] = "version four";
void v4(int dim, pixel *src, pixel *dst)
{
    int i, j, a;
    for (i = 0; i < dim; i++){
	a = dim-1-i;	
	for (j=0; j < dim; j++){
	    dst[RIDX(j, a, dim)] = src[RIDX(i, j, dim)];
	}
    }
}
char v5_descr[] = "version five";
void v5(int dim, pixel *src, pixel *dst)
{
    int i, j, a, b, c, d;
    for (i = 0; i < dim; i+=2){
	a = dim-1-i;
	b = a-1;
	d = i+1;	
	for(j = 0; j < dim; j+=2){
	    c = j+1;
	    dst[RIDX(j, a, dim)] = src[RIDX(i, j, dim)];
	    dst[RIDX(j, b, dim)] = src[RIDX(d, j, dim)];
	    
	    dst[RIDX(c, a, dim)] = src[RIDX(i, c, dim)];
	    dst[RIDX(c, b, dim)] = src[RIDX(d, c, dim)];
	}
    }
}
char v6_descr[] = "version six";
void v6(int dim, pixel *src, pixel *dst)
{
    int i,j,a,b,c,d,e,f,g,h;
    for(i=0; i < dim; i+=4){
	a = dim-1-i;
	b = a-1;
	c = b-1;
	d = c-1;
	f = i+1;
	g = i+2;
	h = i+3;	
	for(j=0; j < dim; j+=4){
	    e = j;

	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    e++;

	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    e++;

	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    e++;

	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	}
    }
}
char v7_descr[] = "version seven";
void v7(int dim, pixel *src, pixel *dst)
{
    int i,j,a,b,c,d,e,f,g,h,k,l,m,p,q,r,s,t;
    for(i=0; i < dim; i+=8){
	a = dim-1-i;
	b = a-1;
	c = b-1;
	d = c-1;
	k = d-1;
	l = k-1;
	m = l-1;
	p = m-1;
	f = i+1;
	g = i+2;
	h = i+3;
	q = i+4;
	r = i+5;
	s = i+6;
	t = i+7;	
	for(j=0; j < dim; j+=8){
	    e = j;

	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
	    e++;

	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
	    e++;

	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
	    e++;

	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
	    e++;
	
	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
	    e++;

	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
	    e++;
	   
	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
	    e++;
	
	    dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	    dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
	}
    }
}
char v8_descr[] = "version eight";
void v8(int dim, pixel *src, pixel *dst)
{
    int i,j,a,b,c,d,e,f,g,h,k,l,m,p,q,r,s,t;
    for(i=0; i < dim; i+=8){
	a = dim-1-i;
	b = a-1;
	c = b-1;
	d = c-1;
	k = d-1;
	l = k-1;
	m = l-1;
	p = m-1;
	f = i+1;
	g = i+2;
	h = i+3;
	q = i+4;
	r = i+5;
	s = i+6;
	t = i+7;	
	for(j=0; j < dim; j+=8){
	    for(e=j; e < j+8; e++){
	    	dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    	dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    	dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    	dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	   	dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    	dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    	dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    	dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
	    }
	}
     }
}
char v9_descr[] = "version nine";
void v9(int dim, pixel *src, pixel *dst)
{
    int a,aa,aaa,b,bb,bbb,c,cc,ccc,d,dd,ddd,e,f,ff,g,gg,h,hh,i,j,k,kk,l,ll,m,mm,p,pp,q,qq,r,s,t;
    for(i=0; i < dim; i+=16){
	a = dim-1-i;
	b = a-1;
	c = b-1;
	d = c-1;
	k = d-1;
	l = k-1;
	m = l-1;
	p = m-1;
	aa = p-1;
	bb = aa-1;
	cc = bb-1;
	dd = cc-1;
	aaa = dd-1;
	bbb = aaa-1;
	ccc = bbb-1;
	ddd = ccc-1;
	f = i+1;
	g = i+2;
	h = i+3;
	q = i+4;
	r = i+5;
	s = i+6;
	t = i+7;
	ff = i+8;
	gg = i+9;
	hh = i+10;
	kk = i+11;
	ll = i+12;
	mm = i+13;
	pp = i+14;
	qq = i+15;	
	for(j=0; j < dim; j+=16){
	    for(e=j; e < j+16; e++){
	        dst[RIDX(e, a, dim)] = src[RIDX(i, e, dim)];
	    	dst[RIDX(e, b, dim)] = src[RIDX(f, e, dim)];
	    	dst[RIDX(e, c, dim)] = src[RIDX(g, e, dim)];
	    	dst[RIDX(e, d, dim)] = src[RIDX(h, e, dim)];
	   	dst[RIDX(e, k, dim)] = src[RIDX(q, e, dim)];
	    	dst[RIDX(e, l, dim)] = src[RIDX(r, e, dim)];
	    	dst[RIDX(e, m, dim)] = src[RIDX(s, e, dim)];
	    	dst[RIDX(e, p, dim)] = src[RIDX(t, e, dim)];
		dst[RIDX(e, aa, dim)] = src[RIDX(ff, e, dim)];
		dst[RIDX(e, bb, dim)] = src[RIDX(gg, e, dim)];
		dst[RIDX(e, cc, dim)] = src[RIDX(hh, e, dim)];
		dst[RIDX(e, dd, dim)] = src[RIDX(kk, e, dim)];
		dst[RIDX(e, aaa, dim)] = src[RIDX(ll, e, dim)];
		dst[RIDX(e, bbb, dim)] = src[RIDX(mm, e, dim)];
		dst[RIDX(e, ccc, dim)] = src[RIDX(pp, e, dim)];
		dst[RIDX(e, ddd, dim)] = src[RIDX(qq, e, dim)];
	    }
	}
    }
}
char v10_descr[] = "version ten"; //get rid of function calls also loop unrolling
void v10(int dim, pixel *src, pixel *dst)// i*dim +e
{
    int i,j,a,e,c,b = dim;
    for(i=0 ; i < b ; i+=16){
	a = b-1-i;
	for(j=0; j < b; j+=16){
	    for(e=j; e < j+16; e+=2){
		c = e*dim+a;		
		dst[c] = src[i*dim+e];
		dst[c-1] = src[(i+1)*dim+e];
		dst[c-2] = src[(i+2)*dim+e];
		dst[c-3] = src[(i+3)*dim+e];
		dst[c-4] = src[(i+4)*dim+e];
		dst[c-5] = src[(i+5)*dim+e];
		dst[c-6] = src[(i+6)*dim+e];
		dst[c-7] = src[(i+7)*dim+e];
		dst[c-8] = src[(i+8)*dim+e];
		dst[c-9] = src[(i+9)*dim+e];
		dst[c-10] = src[(i+10)*dim+e];
		dst[c-11] = src[(i+11)*dim+e];
		dst[c-12] = src[(i+12)*dim+e];
		dst[c-13] = src[(i+13)*dim+e];
		dst[c-14] = src[(i+14)*dim+e];
		dst[c-15] = src[(i+15)*dim+e];

		c = (e+1)*dim+a;		
		dst[c] = src[i*dim+e+1];
		dst[c-1] = src[(i+1)*dim+e+1];
		dst[c-2] = src[(i+2)*dim+e+1];
		dst[c-3] = src[(i+3)*dim+e+1];
		dst[c-4] = src[(i+4)*dim+e+1];
		dst[c-5] = src[(i+5)*dim+e+1];
		dst[c-6] = src[(i+6)*dim+e+1];
		dst[c-7] = src[(i+7)*dim+e+1];
		dst[c-8] = src[(i+8)*dim+e+1];
		dst[c-9] = src[(i+9)*dim+e+1];
		dst[c-10] = src[(i+10)*dim+e+1];
		dst[c-11] = src[(i+11)*dim+e+1];
		dst[c-12] = src[(i+12)*dim+e+1];
		dst[c-13] = src[(i+13)*dim+e+1];
		dst[c-14] = src[(i+14)*dim+e+1];
		dst[c-15] = src[(i+15)*dim+e+1];
	    }
	}
    }
}
char v11_descr[] = "version eleven";
void v11(int dim, pixel *src, pixel *dst)
{
    int i,j,a;
    pixel sauce;
    for(i=0 ; i < dim ; i++){
	a = dim-1-i;
	for(j=0; j < dim; j++){
	    sauce.red = src[RIDX(i,j,dim)].red;
	    sauce.green = src[RIDX(i,j,dim)].green;
	    sauce.blue = src[RIDX(i,j,dim)].blue;
            dst[RIDX(j,a,dim)].red = sauce.red;
	    dst[RIDX(j,a,dim)].green = sauce.green;
    	    dst[RIDX(j,a,dim)].blue = sauce.blue;
	}
    }
}
char v12_descr[] = "version twelve";
void v12(int dim, pixel *src, pixel *dst)
{
    int i,j,a;
    pixel sauce;
    for(i=0 ; i < dim ; i+=8){
	a = dim-1-i;
	for(j=0; j < dim; j+=2){

	    sauce.red = src[RIDX(i,j,dim)].red;
	    sauce.green = src[RIDX(i,j,dim)].green;
	    sauce.blue = src[RIDX(i,j,dim)].blue;
            dst[RIDX(j,a,dim)].red = sauce.red;
	    dst[RIDX(j,a,dim)].green = sauce.green;
    	    dst[RIDX(j,a,dim)].blue = sauce.blue;

	    sauce.red = src[RIDX(i+1,j,dim)].red;
	    sauce.green = src[RIDX(i+1,j,dim)].green;
	    sauce.blue = src[RIDX(i+1,j,dim)].blue;
            dst[RIDX(j,a-1,dim)].red = sauce.red;
	    dst[RIDX(j,a-1,dim)].green = sauce.green;
    	    dst[RIDX(j,a-1,dim)].blue = sauce.blue;

 	    sauce.red = src[RIDX(i+2,j,dim)].red;
	    sauce.green = src[RIDX(i+2,j,dim)].green;
	    sauce.blue = src[RIDX(i+2,j,dim)].blue;
            dst[RIDX(j,a-2,dim)].red = sauce.red;
	    dst[RIDX(j,a-2,dim)].green = sauce.green;
    	    dst[RIDX(j,a-2,dim)].blue = sauce.blue;

	    sauce.red = src[RIDX(i+3,j,dim)].red;
	    sauce.green = src[RIDX(i+3,j,dim)].green;
	    sauce.blue = src[RIDX(i+3,j,dim)].blue;
            dst[RIDX(j,a-3,dim)].red = sauce.red;
	    dst[RIDX(j,a-3,dim)].green = sauce.green;
    	    dst[RIDX(j,a-3,dim)].blue = sauce.blue;
	}
    }
}

/* 
 * naive_flip - The naive baseline version of flip 
 */
char naive_flip_descr[] = "naive_flip: Naive baseline implementation";
void naive_flip(int dim, pixel *src, pixel *dst) 
{
    int i, j;
    for (i = 0; i < dim; i++){        
	for (j = 0; j < dim; j++){
            dst[RIDX_F(i, j, dim)].red   = src[RIDX(i, j, dim)].red;
            dst[RIDX_F(i, j, dim)].green = src[RIDX(i, j, dim)].green;
            dst[RIDX_F(i, j, dim)].blue  = src[RIDX(i, j, dim)].blue;
        }
    }
}

/* 
 * flip - Your current working version of flip
 * IMPORTANT: This is the version you will be graded on
 */
char flip_descr[] = "flip: Current working version";
void flip(int dim, pixel *src, pixel *dst) 
{
    v10(dim, src, dst);
}

/*********************************************************************
 * register_flip_functions - Register all of your different versions
 *     of the flip kernel with the driver by calling the
 *     add_flip_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_flip_functions() 
{
    add_flip_function(&flip, flip_descr);   
    //add_flip_function(&naive_flip, naive_flip_descr);   
    /* ... Register additional test functions here */
    add_flip_function(&version_one, version_one_descr);
    add_flip_function(&v2, v2_descr);
    add_flip_function(&v3, v3_descr);
    add_flip_function(&v4, v4_descr);
    add_flip_function(&v5, v5_descr);
    add_flip_function(&v6, v6_descr);
    add_flip_function(&v7, v7_descr);
    add_flip_function(&v8, v8_descr);
    add_flip_function(&v9, v9_descr);
    add_flip_function(&v10, v10_descr);
    add_flip_function(&v11, v11_descr);
    add_flip_function(&v12, v12_descr);
}


/***************
 * CONVOLVE KERNEL
 **************/
 
/***************************************************************
 * Various typedefs and helper functions for the convolve function
 * You may modify these any way you like.
 **************************************************************/

/* A struct used to compute a pixel value */
typedef struct {
    float red;
    float green;
    float blue;
    float weight;
} pixel_sum;

/******************************************************
 * Your different versions of the convolve kernel go here
 ******************************************************/
char convolve_v1_descr[] = "convolve version 1";
void convolve_v1(int dim, pixel *src, pixel *dst)
{
    int i, j, ii, jj, curI, curJ;
    pixel_sum ps;
    
    for (i = 0; i < dim; i++){
        for (j = 0; j < dim; j++){
            ps.red    = 0.0;
            ps.green  = 0.0;
            ps.blue   = 0.0;
            ps.weight = 0.0;
            for (ii = -2; ii <= 2; ii++){
                for (jj = -2; jj <= 2; jj++){
                    curJ = j+jj;
                    if(curJ<0 || curJ>=dim){
                        continue;
                    }
                    curI = i+ii;
                    if(curI<0 || curI>=dim){
                        continue;
                    }
                    ps.red   += src[RIDX(curI, curJ, dim)].red *   kernel[ii+2][jj+2];
                    ps.green += src[RIDX(curI, curJ, dim)].green * kernel[ii+2][jj+2];
                    ps.blue  += src[RIDX(curI, curJ, dim)].blue *  kernel[ii+2][jj+2];
                    ps.weight += kernel[ii+2][jj+2];
                }
            }
            dst[RIDX(i,j,dim)].red   = (unsigned short)(ps.red/ps.weight);
            dst[RIDX(i,j,dim)].green = (unsigned short)(ps.green/ps.weight);
            dst[RIDX(i,j,dim)].blue  = (unsigned short)(ps.blue/ps.weight);
        }
    }
}
char convolve_v2_descr[] = "convolve version 2"; //emboss bottom right
void convolve_v2(int dim, pixel *src, pixel *dst)
{
    int i, j, ii, jj, curI, curJ, a, b;
    pixel_sum ps;
//Top left corner    
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 2; ii <= 4; ii++){
	i = -2 + ii;		
	for(jj = 2; jj <= 4; jj++){
	    j = -2 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(0,0,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(0,0,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(0,0,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//pixel i=0 j=1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 2; ii <= 4; ii++){
	i = -2 + ii;		
	for(jj = 1; jj <= 4; jj++){
	    j = -1 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(0,1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(0,1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(0,1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

// pixel i=1 j=0
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 1; ii <= 4; ii++){
	i = -1 + ii;		
	for(jj = 2; jj <= 4; jj++){
	    j = -2 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(1,0,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(1,0,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(1,0,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//pixel i=1 j=1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 1; ii <= 4; ii++){
	i = -1 + ii;		
	for(jj = 1; jj <= 4; jj++){
	    j = -1 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(1,1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(1,1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(1,1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//Top right corner
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 2; ii <= 4; ii++){
	i = -2 + ii;		
	for(jj = 0; jj <= 2; jj++){
	    j = dim-3 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(0,dim-1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(0,dim-1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(0,dim-1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//corner i=0 j=94
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 2; ii <= 4; ii++){
	i = -2 + ii;		
	for(jj = 0; jj <= 3; jj++){
	    j = dim-4 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(0,dim-2,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(0,dim-2,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(0,dim-2,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//corner i=1 j=dim-1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 1; ii <= 4; ii++){
	i = -1 + ii;		
	for(jj = 0; jj <= 2; jj++){
	    j = dim-3 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(1,dim-1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(1,dim-1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(1,dim-1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//corner i=1 j=dim-2
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 1; ii <= 4; ii++){
	i = -1 + ii;		
	for(jj = 0; jj <= 3; jj++){
	    j = dim-4 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(1,dim-2,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(1,dim-2,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(1,dim-2,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//pixel i=dim-2 j=0
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 0; ii <= 3; ii++){
	i = dim-4 + ii;		
	for(jj = 2; jj <= 4; jj++){
	    j = -2 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(dim-2,0,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-2,0,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-2,0,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//pixel i=dim-2 j=1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 0; ii <= 3; ii++){
	i = dim-4 + ii;		
	for(jj = 1; jj <= 4; jj++){
	    j = -1 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(dim-2,1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-2,1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-2,1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//pixel i=dim-1 j=1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 0; ii <= 2; ii++){
	i = dim-3 + ii;		
	for(jj = 1; jj <= 4; jj++){
	    j = -1 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(dim-1,1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-1,1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-1,1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

// pixel i=dim-2 j=dim-2
   ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 0; ii <= 3; ii++){
	i = dim-4 + ii;		
	for(jj = 0; jj <= 3; jj++){
	    j = dim-4 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(dim-2,dim-2,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-2,dim-2,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-2,dim-2,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//pixel i=dim-2 j=dim-1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 0; ii <= 3; ii++){
	i = dim-4 + ii;		
	for(jj = 0; jj <= 2; jj++){
	    j = dim-3 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(dim-2,dim-1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-2,dim-1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-2,dim-1,dim)].blue = (unsigned short)(ps.blue/ps.weight);
//pixel i=dim-1 j=dim-2
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 0; ii <= 2; ii++){
	i = dim-3 + ii;		
	for(jj = 0; jj <= 3; jj++){
	    j = dim-4 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(dim-1,dim-2,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-1,dim-2,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-1,dim-2,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//top row
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 0.0;
    	for(ii = 2; ii <= 4; ii++){
	    i = -2 + ii;		
	    for(jj = 0; jj <= 4; jj++){
	    	j = a-2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    	ps.weight += kernel[ii][jj];
	    }
    	}
        dst[RIDX(0,a,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(0,a,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(0,a,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }

//2nd row
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 0.0;
    	for(ii = 1; ii <= 4; ii++){
	    i = -1 + ii;		
	    for(jj = 0; jj <= 4; jj++){
	    	j = a-2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    	ps.weight += kernel[ii][jj];
	    }
    	}
        dst[RIDX(1,a,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(1,a,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(1,a,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }
    

//Bottom left corner
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 0; ii <= 2; ii++){
	i = dim-3 + ii;		
	for(jj = 2; jj <= 4; jj++){
	    j = -2 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(dim-1,0,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-1,0,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-1,0,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//Bottom right corner
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 0.0;
    for(ii = 0; ii <= 2; ii++){
	i = dim-3 + ii;		
	for(jj = 0; jj <= 2; jj++){
	    j = dim-3 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    ps.weight += kernel[ii][jj];
	}
    }
    dst[RIDX(dim-1,dim-1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-1,dim-1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-1,dim-1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//bottom rows
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 0.0;
    	for(ii = 0; ii <= 3; ii++){
	    i = dim-4 + ii;		
	    for(jj = 0; jj <= 4; jj++){
	    	j = a-2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    	ps.weight += kernel[ii][jj];
	    }
    	}
        dst[RIDX(dim-2,a,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(dim-2,a,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(dim-2,a,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 0.0;
    	for(ii = 0; ii <= 2; ii++){
	    i = dim-3 + ii;		
	    for(jj = 0; jj <= 4; jj++){
	    	j = a-2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    	ps.weight += kernel[ii][jj];
	    }
    	}
        dst[RIDX(dim-1,a,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(dim-1,a,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(dim-1,a,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }

//leftmost two rows except for corners
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 0.0;
    	for(ii = 0; ii <= 4; ii++){
	    i = a-2 + ii;		
	    for(jj = 2; jj <= 4; jj++){
	    	j = -2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    	ps.weight += kernel[ii][jj];
	    }
    	}
        dst[RIDX(a,0,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(a,0,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(a,0,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 0.0;
    	for(ii = 0; ii <= 4; ii++){
	    i = a-2 + ii;		
	    for(jj = 1; jj <= 4; jj++){
	    	j = -1 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    	ps.weight += kernel[ii][jj];
	    }
    	}
        dst[RIDX(a,1,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(a,1,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(a,1,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }

//rightmost rows
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 0.0;
    	for(ii = 0; ii <= 4; ii++){
	    i = a-2 + ii;		
	    for(jj = 0; jj <= 3; jj++){
	    	j = dim-4 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    	ps.weight += kernel[ii][jj];
	    }
    	}
        dst[RIDX(a,dim-2,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(a,dim-2,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(a,dim-2,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 0.0;
    	for(ii = 0; ii <= 4; ii++){
	    i = a-2 + ii;		
	    for(jj = 0; jj <= 2; jj++){
	    	j = dim-3 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    	ps.weight += kernel[ii][jj];
	    }
    	}
        dst[RIDX(a,dim-1,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(a,dim-1,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(a,dim-1,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }

    for (i = 2; i < dim-2; i++){
        for (j = 2; j < dim-2; j++){
            ps.red    = 0.0;
            ps.green  = 0.0;
            ps.blue   = 0.0;
            ps.weight = 0.0;
	    for (ii = -2; ii <= 2; ii++){
		curI = i+ii;               
		for (jj = -2; jj <= 2; jj++){
                    curJ = j+jj;
                    ps.red   += src[RIDX(curI, curJ, dim)].red *   kernel[ii+2][jj+2];
                    ps.green += src[RIDX(curI, curJ, dim)].green * kernel[ii+2][jj+2];
                    ps.blue  += src[RIDX(curI, curJ, dim)].blue *  kernel[ii+2][jj+2];
                    ps.weight += kernel[ii+2][jj+2];
                }
            }
            dst[RIDX(i,j,dim)].red   = (unsigned short)(ps.red/ps.weight);
            dst[RIDX(i,j,dim)].green = (unsigned short)(ps.green/ps.weight);
            dst[RIDX(i,j,dim)].blue  = (unsigned short)(ps.blue/ps.weight);
        }
    }
}

char convolve_v3_descr[] = "convolve version 3";
void convolve_v3(int dim, pixel *src, pixel *dst)
{
    int i, j, ii, jj, curI, a;
    pixel_sum ps;

    ps.red = 0.0;
    ps.green = 0.0;
    ps.blue = 0.0;
    for(ii=2 ; ii <= 4; ii++){
	i = -2 + ii;	
 	ps.red += src[RIDX(i, 0, dim)].red * kernel[ii][2];
	ps.green += src[RIDX(i, 0, dim)].green * kernel[ii][2];
	ps.blue += src[RIDX(i, 0, dim)].blue * kernel[ii][2];
	
	ps.red += src[RIDX(i, 1, dim)].red * kernel[ii][3];
	ps.green += src[RIDX(i, 1, dim)].green * kernel[ii][3];
	ps.blue += src[RIDX(i, 1, dim)].blue * kernel[ii][3];

	ps.red += src[RIDX(i, 2, dim)].red * kernel[ii][4];
	ps.green += src[RIDX(i, 2, dim)].green * kernel[ii][4];
	ps.blue += src[RIDX(i, 2, dim)].blue * kernel[ii][4];
    }
//Top left corner    
    
    dst[RIDX(0,0,dim)].red = (unsigned short)(ps.red/-28);
    dst[RIDX(0,0,dim)].green = (unsigned short)(ps.green/-28);
    dst[RIDX(0,0,dim)].blue = (unsigned short)(ps.blue/-28);

//pixel i=0 j=1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = -25.0;
    for(ii = 2; ii <= 4; ii++){
	i = -2 + ii;		
	ps.red += src[RIDX(i, 0, dim)].red * kernel[ii][1];
	ps.green += src[RIDX(i, 0, dim)].green * kernel[ii][1];
	ps.blue += src[RIDX(i, 0, dim)].blue * kernel[ii][1];
	
	ps.red += src[RIDX(i, 1, dim)].red * kernel[ii][2];
	ps.green += src[RIDX(i, 1, dim)].green * kernel[ii][2];
	ps.blue += src[RIDX(i, 1, dim)].blue * kernel[ii][2];

	ps.red += src[RIDX(i, 2, dim)].red * kernel[ii][3];
	ps.green += src[RIDX(i, 2, dim)].green * kernel[ii][3];
	ps.blue += src[RIDX(i, 2, dim)].blue * kernel[ii][3];

	ps.red += src[RIDX(i, 3, dim)].red * kernel[ii][4];
	ps.green += src[RIDX(i, 3, dim)].green * kernel[ii][4];
	ps.blue += src[RIDX(i, 3, dim)].blue * kernel[ii][4];
    }
    dst[RIDX(0,1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(0,1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(0,1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

// pixel i=1 j=0
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = -25.0;
    for(ii = 1; ii <= 4; ii++){
	i = -1 + ii;		
	for(jj = 2; jj <= 4; jj++){
	    j = -2 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(1,0,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(1,0,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(1,0,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//pixel i=1 j=1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = -6.0;
    for(ii = 1; ii <= 4; ii++){
	i = -1 + ii;		
	for(jj = 1; jj <= 4; jj++){
	    j = -1 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(1,1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(1,1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(1,1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//Top right corner
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 1.0;
    for(ii = 2; ii <= 4; ii++){
	i = -2 + ii;		
	for(jj = 0; jj <= 2; jj++){
	    j = dim-3 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(0,dim-1,dim)].red = (unsigned short)(ps.red);
    dst[RIDX(0,dim-1,dim)].green = (unsigned short)(ps.green);
    dst[RIDX(0,dim-1,dim)].blue = (unsigned short)(ps.blue);

//corner i=0 j=dim-1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = -20.0;
    for(ii = 2; ii <= 4; ii++){
	i = -2 + ii;		
	for(jj = 0; jj <= 3; jj++){
	    j = dim-4 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(0,dim-2,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(0,dim-2,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(0,dim-2,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//corner i=1 j=dim-1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 22.0;
    for(ii = 1; ii <= 4; ii++){
	i = -1 + ii;		
	for(jj = 0; jj <= 2; jj++){
	    j = dim-3 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(1,dim-1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(1,dim-1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(1,dim-1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//corner i=1 j=dim-2
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 1.0;
    for(ii = 1; ii <= 4; ii++){
	i = -1 + ii;		
	for(jj = 0; jj <= 3; jj++){
	    j = dim-4 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(1,dim-2,dim)].red = (unsigned short)(ps.red);
    dst[RIDX(1,dim-2,dim)].green = (unsigned short)(ps.green);
    dst[RIDX(1,dim-2,dim)].blue = (unsigned short)(ps.blue);

//pixel i=dim-2 j=0
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = -20.0;
    for(ii = 0; ii <= 3; ii++){
	i = dim-4 + ii;		
	for(jj = 2; jj <= 4; jj++){
	    j = -2 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(dim-2,0,dim)].red = (unsigned short)(ps.red/-20.0);
    dst[RIDX(dim-2,0,dim)].green = (unsigned short)(ps.green/-20.0);
    dst[RIDX(dim-2,0,dim)].blue = (unsigned short)(ps.blue/-20.0);

//pixel i=dim-2 j=1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 1.0;
    for(ii = 0; ii <= 3; ii++){
	i = dim-4 + ii;		
	for(jj = 1; jj <= 4; jj++){
	    j = -1 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(dim-2,1,dim)].red = (unsigned short)(ps.red);
    dst[RIDX(dim-2,1,dim)].green = (unsigned short)(ps.green);
    dst[RIDX(dim-2,1,dim)].blue = (unsigned short)(ps.blue);

//pixel i=dim-1 j=1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 22.0;
    for(ii = 0; ii <= 2; ii++){
	i = dim-3 + ii;		
	for(jj = 1; jj <= 4; jj++){
	    j = -1 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(dim-1,1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-1,1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-1,1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

// pixel i=dim-2 j=dim-2
   ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 8.0;
    for(ii = 0; ii <= 3; ii++){
	i = dim-4 + ii;		
	for(jj = 0; jj <= 3; jj++){
	    j = dim-4 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(dim-2,dim-2,dim)].red = (unsigned short)(ps.red / 8.0);
    dst[RIDX(dim-2,dim-2,dim)].green = (unsigned short)(ps.green / 8.0);
    dst[RIDX(dim-2,dim-2,dim)].blue = (unsigned short)(ps.blue / 8.0);

//pixel i=dim-2 j=dim-1
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 27.0;
    for(ii = 0; ii <= 3; ii++){
	i = dim-4 + ii;		
	for(jj = 0; jj <= 2; jj++){
	    j = dim-3 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(dim-2,dim-1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-2,dim-1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-2,dim-1,dim)].blue = (unsigned short)(ps.blue/ps.weight);
//pixel i=dim-1 j=dim-2
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 27.0;
    for(ii = 0; ii <= 2; ii++){
	i = dim-3 + ii;		
	for(jj = 0; jj <= 3; jj++){
	    j = dim-4 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(dim-1,dim-2,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-1,dim-2,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-1,dim-2,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//top row
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	for(ii = 2; ii <= 4; ii++){
	    i = -2 + ii;		
	    for(jj = 0; jj <= 4; jj++){
	    	j = a-2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    }
    	}
        dst[RIDX(0,a,dim)].red = (unsigned short)(ps.red/-23);
        dst[RIDX(0,a,dim)].green = (unsigned short)(ps.green/-23);
        dst[RIDX(0,a,dim)].blue = (unsigned short)(ps.blue/-23);
    }

//2nd row
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = -3.0;
    	for(ii = 1; ii <= 4; ii++){
	    i = -1 + ii;		
	    for(jj = 0; jj <= 4; jj++){
	    	j = a-2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    }
    	}
        dst[RIDX(1,a,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(1,a,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(1,a,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }
    

//Bottom left corner
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 1.0;
    for(ii = 0; ii <= 2; ii++){
	i = dim-3 + ii;		
	for(jj = 2; jj <= 4; jj++){
	    j = -2 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(dim-1,0,dim)].red = (unsigned short)(ps.red);
    dst[RIDX(dim-1,0,dim)].green = (unsigned short)(ps.green);
    dst[RIDX(dim-1,0,dim)].blue = (unsigned short)(ps.blue);

//Bottom right corner
    ps.red    = 0.0;
    ps.green  = 0.0;
    ps.blue   = 0.0;
    ps.weight = 30.0;
    for(ii = 0; ii <= 2; ii++){
	i = dim-3 + ii;		
	for(jj = 0; jj <= 2; jj++){
	    j = dim-3 + jj;
	    ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	}
    }
    dst[RIDX(dim-1,dim-1,dim)].red = (unsigned short)(ps.red/ps.weight);
    dst[RIDX(dim-1,dim-1,dim)].green = (unsigned short)(ps.green/ps.weight);
    dst[RIDX(dim-1,dim-1,dim)].blue = (unsigned short)(ps.blue/ps.weight);

//bottom rows
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 5.0;
    	for(ii = 0; ii <= 3; ii++){
	    i = dim-4 + ii;		
	    for(jj = 0; jj <= 4; jj++){
	    	j = a-2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    }
    	}
        dst[RIDX(dim-2,a,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(dim-2,a,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(dim-2,a,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 25.0;
    	for(ii = 0; ii <= 2; ii++){
	    i = dim-3 + ii;		
	    for(jj = 0; jj <= 4; jj++){
	    	j = a-2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    }
    	}
        dst[RIDX(dim-1,a,dim)].red = (unsigned short)(ps.red/25.0);
        dst[RIDX(dim-1,a,dim)].green = (unsigned short)(ps.green/25.0);
        dst[RIDX(dim-1,a,dim)].blue = (unsigned short)(ps.blue/25.0);
    }

//leftmost two rows except for corners
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = -23.0;
    	for(ii = 0; ii <= 4; ii++){
	    i = a-2 + ii;		
	    for(jj = 2; jj <= 4; jj++){
	    	j = -2 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    }
    	}
        dst[RIDX(a,0,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(a,0,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(a,0,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = -3.0;
    	for(ii = 0; ii <= 4; ii++){
	    i = a-2 + ii;		
	    for(jj = 1; jj <= 4; jj++){
	    	j = -1 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    }
    	}
        dst[RIDX(a,1,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(a,1,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(a,1,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }

//rightmost rows
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 5.0;
    	for(ii = 0; ii <= 4; ii++){
	    i = a-2 + ii;		
	    for(jj = 0; jj <= 3; jj++){
	    	j = dim-4 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    }
    	}
        dst[RIDX(a,dim-2,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(a,dim-2,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(a,dim-2,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }
    for(a=2; a < dim-2; a++){
        ps.red    = 0.0;
    	ps.green  = 0.0;
   	ps.blue   = 0.0;
    	ps.weight = 25.0;
    	for(ii = 0; ii <= 4; ii++){
	    i = a-2 + ii;		
	    for(jj = 0; jj <= 2; jj++){
	    	j = dim-3 + jj;
	    	ps.red += src[RIDX(i, j, dim)].red * kernel[ii][jj];
	    	ps.green += src[RIDX(i, j, dim)].green * kernel[ii][jj];
	    	ps.blue += src[RIDX(i, j, dim)].blue * kernel[ii][jj];
	    }
    	}
        dst[RIDX(a,dim-1,dim)].red = (unsigned short)(ps.red/ps.weight);
        dst[RIDX(a,dim-1,dim)].green = (unsigned short)(ps.green/ps.weight);
        dst[RIDX(a,dim-1,dim)].blue = (unsigned short)(ps.blue/ps.weight);
    }

    for (i = 2; i < dim-2; i++){
        for (j = 0; j < dim-4; j++){
            ps.red    = 0.0;
            ps.green  = 0.0;
            ps.blue   = 0.0;
	    for (ii = -2; ii <= 2; ii++){
		curI = i+ii;               
		ps.red += src[RIDX(curI, j, dim)].red * kernel[ii+2][0];
		ps.green += src[RIDX(curI, j, dim)].green * kernel[ii+2][0];
		ps.blue += src[RIDX(curI, j, dim)].blue * kernel[ii+2][0];		

		ps.red += src[RIDX(curI, j+1, dim)].red * kernel[ii+2][1];
		ps.green += src[RIDX(curI,j+1, dim)].green * kernel[ii+2][1];
		ps.blue += src[RIDX(curI, j+1, dim)].blue * kernel[ii+2][1];
	
		ps.red += src[RIDX(curI, j+2, dim)].red * kernel[ii+2][2];
		ps.green += src[RIDX(curI, j+2, dim)].green * kernel[ii+2][2];
		ps.blue += src[RIDX(curI, j+2, dim)].blue * kernel[ii+2][2];

		ps.red += src[RIDX(curI, j+3, dim)].red * kernel[ii+2][3];
		ps.green += src[RIDX(curI, j+3, dim)].green * kernel[ii+2][3];
		ps.blue += src[RIDX(curI, j+3, dim)].blue * kernel[ii+2][3];
	
		ps.red += src[RIDX(curI, j+4, dim)].red * kernel[ii+2][4];
		ps.green += src[RIDX(curI, j+4, dim)].green * kernel[ii+2][4];
		ps.blue += src[RIDX(curI, j+4, dim)].blue * kernel[ii+2][4];
            }
            dst[RIDX(i,j+2,dim)].red   = (unsigned short)(ps.red);
            dst[RIDX(i,j+2,dim)].green = (unsigned short)(ps.green);
            dst[RIDX(i,j+2,dim)].blue  = (unsigned short)(ps.blue);
        }
    }
}
/*
 * naive_convolve - The naive baseline version of convolve 
 */
char naive_convolve_descr[] = "naive_convolve: Naive baseline implementation";
void naive_convolve(int dim, pixel *src, pixel *dst) 
{
    int i, j, ii, jj, curI, curJ;
    pixel_sum ps;
    
    for (j = 0; j < dim; j++){
        for (i = 0; i < dim; i++){
            ps.red    = 0.0;
            ps.green  = 0.0;
            ps.blue   = 0.0;
            ps.weight = 0.0;
            for (jj = -2; jj <= 2; jj++){
                for (ii = -2; ii <= 2; ii++){
                    curJ = j+jj;
                    if(curJ<0 || curJ>=dim){
                        continue;
                    }
                    curI = i+ii;
                    if(curI<0 || curI>=dim){
                        continue;
                    }
                    ps.red   += src[RIDX(curI, curJ, dim)].red *   kernel[ii+2][jj+2];
                    ps.green += src[RIDX(curI, curJ, dim)].green * kernel[ii+2][jj+2];
                    ps.blue  += src[RIDX(curI, curJ, dim)].blue *  kernel[ii+2][jj+2];
                    ps.weight += kernel[ii+2][jj+2];
                }
            }
            dst[RIDX(i,j,dim)].red   = (unsigned short)(ps.red/ps.weight);
            dst[RIDX(i,j,dim)].green = (unsigned short)(ps.green/ps.weight);
            dst[RIDX(i,j,dim)].blue  = (unsigned short)(ps.blue/ps.weight);
        }
    }
}

/*
 * convolve - Your current working version of convolve. 
 * IMPORTANT: This is the version you will be graded on
 */
char convolve_descr[] = "convolve: Current working version";
void convolve(int dim, pixel *src, pixel *dst) 
{
    convolve_v2(dim, src, dst);
}

/********************************************************************* 
 * register_convolve_functions - Register all of your different versions
 *     of the convolve kernel with the driver by calling the
 *     add_convolve_function() for each test function.  When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_convolve_functions() {
    add_convolve_function(&convolve, convolve_descr);
    //add_convolve_function(&naive_convolve, naive_convolve_descr);
    /* ... Register additional test functions here */
    add_convolve_function(&convolve_v1, convolve_v1_descr);
    add_convolve_function(&convolve_v2, convolve_v2_descr);
    add_convolve_function(&convolve_v3, convolve_v3_descr);
}

